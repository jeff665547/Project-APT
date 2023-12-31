/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the files COPYING and Copyright.html.  COPYING can be found at the root   *
 * of the source code distribution tree; Copyright.html can be found at the  *
 * root level of an installed copy of the electronic HDF5 document set and   *
 * is linked from the top-level documents page.  It can also be found at     *
 * http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have          *
 * access to either file, you may request a copy from help@hdfgroup.org.     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/***********************************************************
*
* Test program:	 tvlstr
*
* Test the variable length string functionality
*
*************************************************************/

#ifdef OLD_HEADER_FILENAME
#include <iostream.h>
#else
#include <iostream>
#endif
#include <string>

#ifndef H5_NO_NAMESPACE
#ifndef H5_NO_STD
    using std::cerr;
    using std::endl;
#endif  // H5_NO_STD
#endif

#include "testhdf5.h"   // C test header file
#include "H5Cpp.h"      // C++ API header file

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

#include "h5cpputil.h"  // C++ utilility header file

const H5std_string      DATAFILE("tvlstr.h5");
const H5std_string      DATAFILE2("tvlstr2.h5");

// 1-D dataset with fixed dimensions
const int SPACE1_RANK = 1;
const hsize_t SPACE1_DIM1 = 4;

const H5std_string      VLSTR_TYPE("vl_string_type");

// Definitions for the VL re-writing test
const int REWRITE_NDATASETS = 32;

/***********************************************************
*
* Test program:	 tvlstr
*
* Test the Variable-Length String functionality
*
*************************************************************/

void *test_vlstr_alloc_custom(size_t size, void *info);
void test_vlstr_free_custom(void *mem, void *info);

/****************************************************************
**
**  test_vlstr_alloc_custom(): Test VL datatype custom memory
**      allocation routines.  This routine just uses malloc to
**      allocate the memory and increments the amount of memory
**      allocated.
**
****************************************************************/
void *test_vlstr_alloc_custom(size_t size, void *info)
{
    void *ret_value=NULL;   	// Pointer to return
    int *mem_used=(int *)info;  // Get the pointer to the memory used
    size_t extra;           	// Extra space needed

    /*
     *  This weird contortion is required on the DEC Alpha to keep the
     *  alignment correct - QAK
     */
    
    extra=MAX(sizeof(void *),sizeof(size_t));

    if((ret_value=HDmalloc(extra+size))!=NULL) {
        *(size_t *)ret_value=size;
        *mem_used+=size;
    } // end if
    ret_value=((unsigned char *)ret_value)+extra;
    return(ret_value);
}

/****************************************************************
**
**  test_vlstr_free_custom(): Test VL datatype custom memory
**      allocation routines.  This routine just uses free to
**      release the memory and decrements the amount of memory
**      allocated.
**
****************************************************************/
void test_vlstr_free_custom(void *_mem, void *info)
{
    unsigned char *mem;
    int *mem_used=(int *)info;  // Get the pointer to the memory used
    size_t extra;           	// Extra space needed

    /*
     *  This weird contortion is required on the DEC Alpha to keep the
     *  alignment correct - QAK
     */
    
    extra=MAX(sizeof(void *),sizeof(size_t));

    if(_mem!=NULL) {
        mem=((unsigned char *)_mem)-extra;
        *mem_used-=*(size_t *)mem;
        HDfree(mem);
    } // end if
}

/****************************************************************
**
**  test_vlstrings_basic(): Test basic VL string code.
**      Tests simple VL string I/O
**
****************************************************************/
static void
test_vlstrings_basic(void)
{
    const char *wdata[SPACE1_DIM1]= {
        "Four score and seven years ago our forefathers brought forth on this continent a new nation,",
        "conceived in liberty and dedicated to the proposition that all men are created equal.",
        "Now we are engaged in a great civil war,",
        "testing whether that nation or any nation so conceived and so dedicated can long endure."
        };   // Information to write

    // Output message about test being performed
    MESSAGE(5, ("Testing Basic VL String Functionality\n"));

    H5File* file1 = NULL;
    try {
        // Create file.
        file1 = new H5File (DATAFILE, H5F_ACC_TRUNC);

	// Create dataspace for datasets.
	hsize_t	dims1[] = {SPACE1_DIM1};
	DataSpace sid1(SPACE1_RANK, dims1);

	// Create a datatype to refer to.
	StrType tid1(0, H5T_VARIABLE);

	// Create a dataset.
	DataSet dataset(file1->createDataSet("Dataset1", tid1, sid1));

	// Write dataset to disk.
	dataset.write(wdata, tid1);

	// Create H5S_SCALAR data space.
	DataSpace dataspace;

	DataSet dataset2(file1->createDataSet("Dataset2", tid1, dataspace));

	char *wdata2 = (char*)HDcalloc(65534, sizeof(char));
	HDmemset(wdata2, 'A', 65533);

	dataset2.write(&wdata2, tid1);

	dataspace.close();
	dataset2.close();
	HDfree(wdata2);

	// Change to the custom memory allocation routines for reading VL string.
	DSetMemXferPropList xfer;
	int mem_used=0;	// Memory used during allocation
	xfer.setVlenMemManager(test_vlstr_alloc_custom, &mem_used, test_vlstr_free_custom, &mem_used);

	// Make certain the correct amount of memory will be used.
	hsize_t vlsize = dataset.getVlenBufSize(tid1, sid1);

	// Count the actual number of bytes used by the strings.
	int     str_used;   // String data in memory
	hsize_t i;	// counting variable
	for (i=0,str_used=0; i<SPACE1_DIM1; i++)
	    str_used+=HDstrlen(wdata[i])+1;

	// Compare against the strings actually written.
	verify_val((int)vlsize,str_used,"DataSet::getVlenBufSize", __LINE__, __FILE__);

	// Read dataset from disk.
	char *rdata[SPACE1_DIM1];   // Information read in
	dataset.read(rdata, tid1, DataSpace::ALL, DataSpace::ALL, xfer);

	// Make certain the correct amount of memory has been used.
	VERIFY(mem_used, str_used, "H5Dread");

	// Compare data read in.
	for (i=0; i<SPACE1_DIM1; i++) {
	    if(HDstrlen(wdata[i])!=strlen(rdata[i])) {
		TestErrPrintf("VL data length don't match!, strlen(wdata[%d])=%d, strlen(rdata[%d])=%d\n",(int)i,(int)strlen(wdata[i]),(int)i,(int)strlen(rdata[i]));
		continue;
	    } // end if
	    if( HDstrcmp(wdata[i],rdata[i]) != 0 ) {
		TestErrPrintf("VL data values don't match!, wdata[%d]=%s, rdata[%d]=%s\n",(int)i,wdata[i],(int)i,rdata[i]);
		continue;
	    } // end if
	} // end for

	// Reclaim the read VL data.
	DataSet::vlenReclaim((void *)rdata, tid1, sid1, xfer);

	// Make certain the VL memory has been freed.
	VERIFY(mem_used, 0, "DataSet::vlenReclaim");

	// Close objects and file.
	dataset.close();
	tid1.close();
	sid1.close();
	xfer.close();
	file1->close();
    } // end try

    // Catch all exceptions.
    catch (Exception E)
    {
	issue_fail_msg(E.getCFuncName(), __LINE__, __FILE__, E.getCDetailMsg());
        if (file1 != NULL) // clean up
            delete file1;
    }
} // end test_vlstrings_basic()

/****************************************************************
**
**  test_vlstrings_special(): Test VL string code for special
**      string cases, nil and zero-sized.
**
****************************************************************/
static void
test_vlstrings_special(void)
{
    const char *wdata[SPACE1_DIM1] = {"one", "two", "", "four"};
    const char *wdata2[SPACE1_DIM1] = {NULL, NULL, NULL, NULL};
    char *rdata[SPACE1_DIM1];   // Information read in

    // Output message about test being performed.
    MESSAGE(5, ("Testing Special VL Strings\n"));

    try {
	// Create file.
	H5File file1(DATAFILE, H5F_ACC_TRUNC);

        // Create dataspace for datasets.
        hsize_t dims1[] = {SPACE1_DIM1};
        DataSpace sid1(SPACE1_RANK, dims1);

	// Create a datatype to refer to.
	StrType tid1(0, H5T_VARIABLE);

	// Create a dataset.
	DataSet dataset(file1.createDataSet("Dataset3", tid1, sid1));

	// Read from dataset before writing data.
	dataset.read(rdata, tid1);

	// Check data read in.
	hsize_t i;      	// counting variable
	for (i=0; i<SPACE1_DIM1; i++)
	    if(rdata[i]!=NULL)
		TestErrPrintf("VL doesn't match!, rdata[%d]=%p\n",(int)i,rdata[i]);

	// Write dataset to disk.
	dataset.write(wdata, tid1);

	// Read dataset from disk.
	dataset.read(rdata, tid1);

	// Compare data read in.
	for (i=0; i<SPACE1_DIM1; i++) {
	    if(HDstrlen(wdata[i])!=strlen(rdata[i])) {
		TestErrPrintf("VL data length don't match!, strlen(wdata[%d])=%d, strlen(rdata[%d])=%d\n",(int)i,(int)strlen(wdata[i]),(int)i,(int)strlen(rdata[i]));
		continue;
	    } // end if
	    if( HDstrcmp(wdata[i],rdata[i]) != 0 ) {
		TestErrPrintf("VL data values don't match!, wdata[%d]=%s, rdata[%d]=%s\n",(int)i,wdata[i],(int)i,rdata[i]);
		continue;
	    } // end if
	} // end for

	// Reclaim the read VL data.
	DataSet::vlenReclaim((void *)rdata, tid1, sid1);

	// Close Dataset.
	dataset.close();

	// Create another dataset to test nil strings.
	DSetCreatPropList dcpl;

	// Set the fill value for the second dataset.
	char *fill = NULL;	// Fill value
	dcpl.setFillValue(tid1, &fill);

	dataset = file1.createDataSet("Dataset4", tid1, sid1, dcpl);

	// Close dataset creation property list.
	dcpl.close();

	// Read from dataset before writing data.
	dataset.read(rdata, tid1);

	// Check data read in.
	for (i=0; i<SPACE1_DIM1; i++)
	  if(rdata[i]!=NULL)
	    TestErrPrintf("VL doesn't match!, rdata[%d]=%p\n",(int)i,rdata[i]);

	// Try to write nil strings to disk.
	dataset.write(wdata2, tid1);

	// Read nil strings back from disk.
	dataset.read(rdata, tid1);

	// Check data read in.
	for (i=0; i<SPACE1_DIM1; i++)
	  if(rdata[i]!=NULL)
	    TestErrPrintf("VL doesn't match!, rdata[%d]=%p\n",(int)i,rdata[i]);

	// Close objects and file.
	dataset.close();
	tid1.close();
	sid1.close();
	file1.close();
    } // end try

    // Catch all exceptions.
    catch (Exception E)
    {
	issue_fail_msg(E.getCFuncName(), __LINE__, __FILE__, E.getCDetailMsg());
    }
} // test_vlstrings_special

/****************************************************************
**
**  test_vlstring_type(): Test VL string type.
**      Tests if VL string is treated as string.
**
****************************************************************/
static void test_vlstring_type(void)
{
    // Output message about test being performed.
    MESSAGE(5, ("Testing VL String type\n"));

    H5File* file1 = NULL;
    try {
	// Open file.
	file1 = new H5File(DATAFILE, H5F_ACC_RDWR);

	// Create a datatype to refer to.
	StrType vlstr_type(PredType::C_S1);

	// Change padding and verify it.
	vlstr_type.setStrpad(H5T_STR_NULLPAD);
	H5T_str_t pad = vlstr_type.getStrpad();
	verify_val(pad, H5T_STR_NULLPAD, "StrType::getStrpad", __LINE__, __FILE__);

	// Convert to variable-length string.
	vlstr_type.setSize(H5T_VARIABLE);

	// Check if datatype is VL string.
	H5T_class_t type_class = vlstr_type.getClass();
	verify_val(type_class, H5T_STRING, "DataType::getClass", __LINE__, __FILE__);
	bool is_variable_str = vlstr_type.isVariableStr();
	verify_val(is_variable_str, true, "DataType::isVariableStr", __LINE__, __FILE__);

	// Check default character set and padding.
	H5T_cset_t cset = vlstr_type.getCset();
	verify_val(cset, H5T_CSET_ASCII, "StrType::getCset", __LINE__, __FILE__);
	pad = vlstr_type.getStrpad();
	verify_val(pad, H5T_STR_NULLPAD, "StrType::getStrpad", __LINE__, __FILE__);

	// Commit variable-length string datatype to storage.
	vlstr_type.commit(*file1, VLSTR_TYPE);

	// Close datatype.
	vlstr_type.close();

	// Try opening datatype again.
	vlstr_type = file1->openStrType(VLSTR_TYPE);

	// Close datatype and file.
	vlstr_type.close();
	file1->close();

	// Open file.
	file1 = new H5File(DATAFILE, H5F_ACC_RDWR);

    //fid = H5Fopen(DATAFILE.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

	// Open the variable-length string datatype just created
	vlstr_type.setId((file1->openStrType(VLSTR_TYPE)).getId());
    //tid_vlstr = H5Topen(fid, VLSTR_TYPE.c_str());

	// Verify character set and padding
	cset = vlstr_type.getCset();
	verify_val(cset, H5T_CSET_ASCII, "StrType::getCset", __LINE__, __FILE__);
    //cset = H5Tget_cset(tid_vlstr);
	pad = vlstr_type.getStrpad();
    //pad = H5Tget_strpad(tid_vlstr);
	verify_val(pad, H5T_STR_NULLPAD, "StrType::getStrpad", __LINE__, __FILE__);

	// Close datatype and file
	vlstr_type.close();
	file1->close();
    } // end try

    // Catch all exceptions.
    catch (Exception E)
    {
        issue_fail_msg(E.getCFuncName(), __LINE__, __FILE__, E.getCDetailMsg());
    }
} // end test_vlstring_type()

/****************************************************************
**
**  test_compact_vlstring(): Test code for storing VL strings in
**      compact datasets.
**
****************************************************************/
static void
test_compact_vlstring(void)
{
    const char *wdata[SPACE1_DIM1] = {"one", "two", "three", "four"};
    char *rdata[SPACE1_DIM1];   // Information read in

	// Output message about test being performed
    MESSAGE(5, ("Testing VL Strings in compact dataset\n"));

    try {
	// Create file
	H5File file1(DATAFILE, H5F_ACC_TRUNC);

	// Create dataspace for datasets
        hsize_t dims1[] = {SPACE1_DIM1};
        DataSpace sid1(SPACE1_RANK, dims1);

	// Create a datatype to refer to
	StrType tid1(0, H5T_VARIABLE);

	// Create dataset create property list and set layout
	DSetCreatPropList plist;
	plist.setLayout(H5D_COMPACT);

	// Create a dataset
	DataSet dataset(file1.createDataSet("Dataset5", tid1, sid1, plist));

	// Write dataset to disk
	dataset.write(wdata, tid1);

	// Read dataset from disk
	dataset.read(rdata, tid1);

	// Compare data read in
	hsize_t i;
	for (i=0; i<SPACE1_DIM1; i++) {
	    if (HDstrlen(wdata[i])!=strlen(rdata[i])) {
		TestErrPrintf("VL data length don't match!, strlen(wdata[%d])=%d, strlen(rdata[%d])=%d\n",(int)i,(int)strlen(wdata[i]),(int)i,(int)strlen(rdata[i]));
		continue;
	    } // end if
	    if (HDstrcmp(wdata[i],rdata[i]) != 0) {
		TestErrPrintf("VL data values don't match!, wdata[%d]=%s, rdata[%d]=%s\n",(int)i,wdata[i],(int)i,rdata[i]);
		continue;
	    } // end if
	} // end for

	// Reclaim the read VL data
	DataSet::vlenReclaim((void *)rdata, tid1, sid1);

	// Close objects and file
	dataset.close();
	tid1.close();
	sid1.close();
	plist.close();
	file1.close();
    } // end try

    // Catch all exceptions.
    catch (Exception E)
    {
        issue_fail_msg(E.getCFuncName(), __LINE__, __FILE__, E.getCDetailMsg());
    }
}   // test_compact_vlstrings

/****************************************************************
**
**  test_write_vl_string_attribute(): Test basic VL string code.
**      Tests writing VL strings as attributes
**
****************************************************************/
// String for testing attributes
static char *string_att_write=NULL;

// Info for a string attribute
const H5std_string ATTRSTR_NAME("String_attr");
const H5std_string ATTRSTR_DATA("String Attribute");

static void 
test_write_vl_string_attribute(void)
{
    // Output message about test being performed
    MESSAGE(5, ("Testing writing VL String as attributes\n"));

    try {
	// Open the file
	H5File file1(DATAFILE, H5F_ACC_RDWR);

	// Create a datatype to refer to.
	StrType tid1(0, H5T_VARIABLE);

	// Open the root group.
	Group root = file1.openGroup("/");

	// Create dataspace for the attribute.
	DataSpace att_space (H5S_SCALAR);

	// Create an attribute for the root group.
	Attribute gr_attr = root.createAttribute(ATTRSTR_NAME, tid1, att_space);

	// Write data to the attribute.
	gr_attr.write(tid1, ATTRSTR_DATA);

	// Read and verify the attribute string as a string of chars.
	char *string_att_check;
	gr_attr.read(tid1, &string_att_check);
	if(HDstrcmp(string_att_check, ATTRSTR_DATA.c_str())!=0)
	    TestErrPrintf("Line %d: Attribute data different: ATTRSTR_DATA=%s,string_att_check=%s\n",__LINE__, ATTRSTR_DATA.c_str(), string_att_check);

	HDfree(string_att_check);  // note: no need for std::string test

	// Read and verify the attribute string as an std::string.
	H5std_string read_str;
	gr_attr.read(tid1, read_str);
	if (read_str != ATTRSTR_DATA)
	    TestErrPrintf("Line %d: Attribute data different: ATTRSTR_DATA=%s,read_str=%s\n",__LINE__, ATTRSTR_DATA.c_str(), read_str.c_str());

	// Close group's attribute.
	gr_attr.close();

	// Test creating a "large" sized string attribute
	gr_attr = root.createAttribute("test_scalar_large", tid1, att_space);

	string_att_write = (char*)HDcalloc(8192, sizeof(char));
	HDmemset(string_att_write, 'A', 8191);

	// Write data to the attribute.
	gr_attr.write(tid1, &string_att_write);

	gr_attr.read(tid1, &string_att_check);

	if(HDstrcmp(string_att_check,string_att_write)!=0)
	    TestErrPrintf("VL string attributes don't match!, string_att_write=%s, string_att_check=%s\n",string_att_write,string_att_check);

	HDfree(string_att_check);

	gr_attr = root.openAttribute(ATTRSTR_NAME);

	// The attribute string written is freed below, in the 
	// test_read_vl_string_attribute() test

	// Close attribute and file
	gr_attr.close();
	file1.close();
    } // end try block

    // Catch all exceptions.
    catch (Exception E) {
	issue_fail_msg(E.getCFuncName(), __LINE__, __FILE__, E.getCDetailMsg());
    }
}   // test_string_attr()

/****************************************************************
**
**  test_read_vl_string_attribute(): Test basic VL string code.
**      Tests reading VL strings from attributes
**
****************************************************************/
static void test_read_vl_string_attribute(void)
{
    char *string_att_check;

    try {
	// Open file
	H5File file1(DATAFILE, H5F_ACC_RDONLY);

	// Create a datatype to refer to.
	StrType tid1(0, H5T_VARIABLE);

	Group root = file1.openGroup("/");

	Attribute att = root.openAttribute(ATTRSTR_NAME);

	// Test reading "normal" sized string attribute
	att.read(tid1, &string_att_check);

	if(HDstrcmp(string_att_check,ATTRSTR_DATA.c_str())!=0)
	    TestErrPrintf("VL string attributes don't match!, string_att=%s, string_att_check=%s\n",ATTRSTR_DATA.c_str(),string_att_check);

	HDfree(string_att_check);

	// Close this attribute.
	att.close();

	// Test reading "large" sized string attribute
	att = root.openAttribute("test_scalar_large");
	att.read(tid1, &string_att_check);

	if(HDstrcmp(string_att_check,string_att_write)!=0)
	    TestErrPrintf("VL string attributes don't match!, string_att_write=%s, string_att_check=%s\n",string_att_write,string_att_check);

	HDfree(string_att_check);
	HDfree(string_att_write);   // Free string allocated in test_write_vl_string_attribute

	// Close objects and file.
	att.close();
	tid1.close();
	root.close();
	file1.close();
    } // end try

    // Catch all exceptions.
    catch (Exception E) {
	issue_fail_msg(E.getCFuncName(), __LINE__, __FILE__, E.getCDetailMsg());
    }
} // test_read_vl_string_attribute

/* Helper routine for test_vl_rewrite() */
static void write_scalar_dset(H5File& file, DataType& type, DataSpace& space, char *name, char *data)
{
    DataSet dset;
    try {
	dset = file.createDataSet(name, type, space);
	dset.write(&data, type, space, space);
	dset.close();
    } // end try
    catch (FileIException ferr) {
	throw;
    }
    catch (DataSetIException derr) {
	throw;
    }
}

/* Helper routine for test_vl_rewrite() */
static void read_scalar_dset(H5File& file, DataType& type, DataSpace& space, char *name, char *data)
{
    char *data_read;
    DataSet dset;
    try {
	dset = file.openDataSet(name);
	dset.read(&data_read, type, space, space);
	dset.close();

	if(HDstrcmp(data, data_read))
	    TestErrPrintf("Expected %s for dataset %s but read %s\n", data, name, data_read);

	HDfree(data_read);
    } // end try
    catch (FileIException ferr) {
	throw;
    }
    catch (DataSetIException derr) {
	throw;
    }
}

/****************************************************************
**
**  test_vl_rewrite(): Test basic VL string code.
**      Tests I/O on VL strings when lots of objects in the file
**      have been linked/unlinked.
**
****************************************************************/
static void test_vl_rewrite(void)
{
    try {
	// Create the files.
	H5File file1(DATAFILE, H5F_ACC_TRUNC);
	H5File file2(DATAFILE2, H5F_ACC_TRUNC);

	// Create the VL string datatype.
	StrType type(0, H5T_VARIABLE);

	// Create dataspace for the attribute.
	DataSpace space (H5S_SCALAR);

	// Create in file 1.
	int i;
	char name[256]; 	// Buffer for names & data
	for (i=0; i<REWRITE_NDATASETS; i++) {
	    sprintf(name, "/set_%d", i);
	    write_scalar_dset(file1, type, space, name, name);
	}

	// Effectively copy data from file 1 to 2.
	for (i=0; i<REWRITE_NDATASETS; i++) {
	    sprintf(name, "/set_%d", i);
	    read_scalar_dset(file1, type, space, name, name);
	    write_scalar_dset(file2, type, space, name, name);
	}

	// Read back from file 2.
	for (i=0; i<REWRITE_NDATASETS; i++) {
	    sprintf(name, "/set_%d", i);
	    read_scalar_dset(file2, type, space, name, name);
	}

	// Remove from file 2.
	for (i=0; i<REWRITE_NDATASETS; i++) {
	    sprintf(name, "/set_%d", i);
	    file2.unlink(name);
	}

	// Effectively copy from file 1 to file 2.
	for (i=0; i<REWRITE_NDATASETS; i++) {
	    sprintf(name, "/set_%d", i);
	    read_scalar_dset(file1, type, space, name, name);
	    write_scalar_dset(file2, type, space, name, name);
	}

	// Close objects and file.
	type.close();
	space.close();
	file1.close();
	file2.close();
    } // end try

    // Catch all exceptions.
    catch (Exception E) {
	issue_fail_msg(E.getCFuncName(), __LINE__, __FILE__, E.getCDetailMsg());
    }
} // end test_vl_rewrite()

/****************************************************************
**
**  test_vlstrings(): Main VL string testing routine.
**
****************************************************************/
void test_vlstrings(void)
{
    // Output message about test being performed
    MESSAGE(5, ("Testing Variable-Length Strings\n"));

    // These tests use the same file
    // Test basic VL string datatype
    test_vlstrings_basic();
    test_vlstrings_special();
    test_vlstring_type();
    test_compact_vlstring();

    // Test using VL strings in attributes
    test_write_vl_string_attribute();
    test_read_vl_string_attribute();

    // Test writing VL datasets in files with lots of unlinking
    test_vl_rewrite();

}   // test_vlstrings()


/*-------------------------------------------------------------------------
 * Function:	cleanup_vlstrings
 *
 * Purpose:	Cleanup temporary test files
 *
 * Return:	none
 *
 * Programmer:	Quincey Koziol
 *              September 10, 1999
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
void
cleanup_vlstrings(void)
{
    HDremove(DATAFILE.c_str());
    HDremove(DATAFILE2.c_str());
}


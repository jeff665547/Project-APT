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

/*****************************************************************************
   FILE
   ttypes.cpp - HDF5 C++ testing the general data type functionality

 ***************************************************************************/

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

const H5std_string      DATAFILE("ttypes.h5");

#define NTESTS	1

// Number of elements in each test
#define NTESTELEM	100000

// Define if you want to see a count of overflows
#undef SHOW_OVERFLOWS

/*
 * Offset from alinged memory returned by malloc().  This can be used to test
 * that type conversions handle non-aligned buffers correctly.
 */

#define ALIGNMENT	1

/*
 * Define if you want to test alignment code on a machine that doesn't
 * normally require alignment. When set, all native data types must be aligned
 * on a byte boundary equal to the data size.
 */

#define TEST_ALIGNMENT

// Alignment test stuff
#ifdef TEST_ALIGNMENT
#define H5T_PACKAGE
#include "H5Tpkg.h"
#endif
#define SET_ALIGNMENT(TYPE,VAL) \
    H5T_NATIVE_##TYPE##_ALIGN_g=MAX(H5T_NATIVE_##TYPE##_ALIGN_g, VAL)

const char *FILENAME[] = {
    "dtypes1.h5",
    "dtypes2.h5",
    "dtypes3.h5",
    NULL
};

/*
 * Count up or down depending on whether the machine is big endian or little
 * endian.  If local variable `endian' is H5T_ORDER_BE then the result will
 * be I, otherwise the result will be Z-(I+1).
 */

#define ENDIAN(Z,I)	(H5T_ORDER_BE==endian?(I):(Z)-((I)+1))


typedef enum flt_t {
    FLT_FLOAT, FLT_DOUBLE, FLT_LDOUBLE, FLT_OTHER
} flt_t;

typedef enum int_t {
    INT_CHAR, INT_UCHAR, INT_SHORT, INT_USHORT, INT_INT, INT_UINT,
    INT_LONG, INT_ULONG, INT_LLONG, INT_ULLONG, INT_OTHER
} int_t;

// Count the number of overflows
#ifdef SHOW_OVERFLOWS
static int noverflows_g = 0;
#endif

// Skip overflow tests if non-zero

// Don't use hardware conversions if set

// Count opaque conversions

/*
 * Although we check whether a floating point overflow generates a SIGFPE and
 * turn off overflow tests in that case, it might still be possible for an
 * overflow condition to occur.  Once a SIGFPE is raised the program cannot
 * be allowed to continue (cf. Posix signals) so in order to recover from a
 * SIGFPE we run tests that might generate one in a child process.
 */

#if defined(H5_HAVE_FORK) && defined(H5_HAVE_WAITPID)
#   define HANDLE_SIGFPE
#endif

// Allocates memory aligned on a certain boundary.
#define aligned_malloc(Z)	((void*)((char*)malloc(ALIGNMENT+Z)+ALIGNMENT))
#define aligned_free(M)		free((char*)(M)-ALIGNMENT)


/*-------------------------------------------------------------------------
 * Function:	overflow_handler
 *
 * Purpose:	Gets called for all data type conversion overflows.
 *
 * Return:	Success:	0
 *
 *		Failure:	-1
 *
 * Programmer:	Robb Matzke
 *              Tuesday, July  7, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
#ifdef SHOW_OVERFLOWS
static herr_t
overflow_handler(hid_t UNUSED src_id, hid_t UNUSED dst_id,
		 void UNUSED *src_buf, void UNUSED *dst_buf)
{
    noverflows_g++;
    return -1;
}
#endif


/*-------------------------------------------------------------------------
 * Function:    test_classes
 *
 * Purpose:     Test type classes
 *
 * Return:      Success:        0
 *
 *              Failure:        number of errors
 *
 * Programmer:  Binh-Minh Ribler (using C version)
 *		January, 2007
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static void test_classes(void)
{
    SUBTEST("PredType::getClass()");
    try {
	// PredType::NATIVE_INT should be in H5T_INTEGER class
	H5T_class_t tcls = PredType::NATIVE_INT.getClass();
	if (H5T_INTEGER!=tcls) {
	    H5_FAILED();
	    puts("    Invalid type class for H5T_NATIVE_INT");
	}

	// PredType::NATIVE_DOUBLE should be in H5T_FLOAT class
	tcls = PredType::NATIVE_DOUBLE.getClass();
	if (H5T_FLOAT!=tcls) {
	    H5_FAILED();
	    puts("    Invalid type class for H5T_NATIVE_DOUBLE");
	}
	PASSED();
    }   // end of try block
    catch (DataTypeIException E) { 
	issue_fail_msg(E.getCFuncName(), __LINE__, __FILE__, E.getCDetailMsg());
    }
}

/*-------------------------------------------------------------------------
 * Function:    test_copy
 *
 * Purpose:     Test datatype copy functionality
 *
 * Return:      Success:        0
 *
 *              Failure:        number of errors
 *
 * Programmer:  Binh-Minh Ribler (using C version)
 *		January, 2007
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static void test_copy(void)
{
    SUBTEST("DataType::copy() and DataType::operator=");
    try {
	// Test copying from a predefined datatype using DataType::operator=
	DataType assigned_type;
	assigned_type = PredType::NATIVE_SHORT;

        // Test copying a predefined type using DataType::copy
	DataType copied_type;
        copied_type.copy (PredType::STD_B8LE);

	// Test copying a user-defined type using DataType::operator=
	DataType assigned_usertype;
	assigned_usertype = copied_type;

	// Test copying from a user-defined datatype using DataType::copy
	DataType copied_usertype;
	copied_usertype.copy(copied_type);

        // Test copying a user-defined int type using DataType::operator=
        IntType orig_int(PredType::STD_B8LE);
        DataType generic_type;
        generic_type = orig_int;

        // Test copying an integer predefined type
        IntType new_int_type(PredType::STD_B8LE);

        // Test copying an int predefined type using DataType::operator=
        IntType another_int_type;
        another_int_type = new_int_type;

	PASSED();
    } 
    catch (DataTypeIException E) { 
	issue_fail_msg(E.getCFuncName(), __LINE__, __FILE__, E.getCDetailMsg());
    }
}


/*-------------------------------------------------------------------------
 * Function:	test_query
 *
 * Purpose:	Tests query functions of compound and enumeration types.
 *
 * Return:	Success: 	0
 * 	
 *		Failure:	number of errors
 *
 * Programmer:	Binh-Minh Ribler (use C version)
 *		January, 2007
 *  
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */

const H5std_string CompT_NAME("Compound_type");
const H5std_string EnumT_NAME("Enum_type");

static void test_query(void)
{
    typedef struct {
	int    a;
	float  b;
	long   c;
	double d;
    } s_type_t;
    short	enum_val;

    // Output message about test being performed
    SUBTEST("Query functions of compound and enumeration types");
    try
    {
	// Create File
	H5File file(FILENAME[2], H5F_ACC_TRUNC);

	// Create a compound datatype
	CompType tid1(sizeof(s_type_t));

	tid1.insertMember("a", HOFFSET(s_type_t, a), PredType::NATIVE_INT);
	tid1.insertMember("b", HOFFSET(s_type_t, b), PredType::NATIVE_FLOAT);
	tid1.insertMember("c", HOFFSET(s_type_t, c), PredType::NATIVE_LONG);
	tid1.insertMember("d", HOFFSET(s_type_t, d), PredType::NATIVE_DOUBLE);

	// Create a enumerate datatype
	EnumType tid2(sizeof(short));

	tid2.insert("RED", (enum_val=0,&enum_val));
	tid2.insert("GREEN", (enum_val=1,&enum_val));
	tid2.insert("BLUE", (enum_val=2,&enum_val));
	tid2.insert("ORANGE", (enum_val=3,&enum_val));
	tid2.insert("YELLOW", (enum_val=4,&enum_val));

	// Query member number and member index by name, for compound type
	int nmembs = tid1.getNmembers();
	verify_val(nmembs, 4, "CompType::getNmembers()", __LINE__, __FILE__);

	int index = tid1.getMemberIndex("c");
	verify_val(index, 2, "CompType::getMemberIndex()", __LINE__, __FILE__);

	// Query member number and member index by name, for enumeration type.
	nmembs = tid2.getNmembers();
	verify_val(nmembs, 5, "EnumType::getNmembers()", __LINE__, __FILE__);

	index = tid2.getMemberIndex("ORANGE");
	verify_val(index, 3, "EnumType::getMemberIndex()", __LINE__, __FILE__);

	// Commit compound datatype and close it
	tid1.commit(file, CompT_NAME);
	tid1.close();

	// Commit enumeration datatype and close it
	tid2.commit(file, EnumT_NAME);
	tid2.close();

	// Open the datatype for query
	tid1 = file.openCompType(CompT_NAME);

	tid2 = file.openEnumType(EnumT_NAME);

	// Query member number and member index by name, for compound type
	nmembs = tid1.getNmembers();
	verify_val(nmembs, 4, "CompType::getNmembers()", __LINE__, __FILE__);

	index = tid1.getMemberIndex("c");
	verify_val(index, 2, "CompType::getMemberIndex()", __LINE__, __FILE__);

	// Query member number and member index by name, for enumeration type
	nmembs = tid2.getNmembers();
	verify_val(nmembs, 5, "EnumType::getNmembers()", __LINE__, __FILE__);

	index = tid2.getMemberIndex("ORANGE");
	verify_val(index, 3, "EnumType::getMemberIndex()", __LINE__, __FILE__);

	// Close data types and file
	tid1.close();
	tid2.close();
	file.close();

    PASSED();
    }   // end of try block
    catch (Exception E) {
        issue_fail_msg(E.getCFuncName(), __LINE__, __FILE__, E.getCDetailMsg());
    }
}   // test_query
 

/*-------------------------------------------------------------------------
 * Function:	test_transient
 *
 * Purpose:	Tests transient data types.
 *
 * Return:	Success:	0
 *
 *		Failure:	number of errors
 *
 * Programmer:	Binh-Minh Ribler (use C version)
 *		January, 2007
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
const H5std_string filename1("dtypes1.h5");
static void test_transient ()
{
    static hsize_t	ds_size[2] = {10, 20};
    
    SUBTEST("Transient data types");
    try {

	H5File file(filename1, H5F_ACC_TRUNC);
	DataSpace space(2, ds_size, ds_size);

	// Predefined types cannot be modified or closed
	// PredType::NATIVE_INT is a constant and cannot make a call, 
	// don't need these tests

	// Copying a predefined type results in a modifiable copy
	IntType type(PredType::NATIVE_INT);
	type.setPrecision(256);

	// It should not be possible to create an attribute for a transient 
	// type 
	try {
	    Attribute attr(type.createAttribute("attr1", PredType::NATIVE_INT, space));
	    // Should FAIL but didn't, so throw an invalid action exception
	    throw InvalidActionException("H5Object::createAttribute", "Attempted to commit a predefined datatype.");
	} catch (AttributeIException err) {}

	// Create a dataset from a transient data type
	type.copy(PredType::NATIVE_INT);
	DataSet dset(file.createDataSet("dset1", type, space));

	// The type returned from a dataset should not be modifiable
	IntType itype(dset);
	try {
	    itype.setPrecision(256);

	    // Should FAIL but didn't, so throw an invalid action exception
	    throw InvalidActionException("PredType::setPrecision", "Dataset data types should not be modifiable!");
	} catch (DataTypeIException err) {}

	itype.close();

	//
 	// Get the dataset data type by applying H5Tcopy() to the dataset. The
 	// result should be modifiable.
 	///
	itype.copy(dset);
	itype.setPrecision(256);

	//
 	// Close the dataset and reopen it, testing that its type is still
 	// read-only. <--- how come modifiable below?
 	///
	dset.close();
    //if (H5Dclose (dset)<0) printf("goto error in C\n");
	dset = file.openDataSet("dset1");
    //if ((dset=H5Dopen (file, "dset1"))<0) printf("goto error in C\n");

	//
 	// Get the dataset data type by applying H5Tcopy() to the dataset. The
 	// result should be modifiable.
 	///
	itype.copy(dset);
	itype.setPrecision(256);
	itype.close();
    
	// Close objects and file.
	dset.close();
	file.close();
	type.close();
	space.close();
	PASSED();
    }   // end of try block
    catch (Exception E) {
        issue_fail_msg(E.getCFuncName(), __LINE__, __FILE__, E.getCDetailMsg());
    }
}   // test_transient


/*-------------------------------------------------------------------------
 * Function:	test_named
 *
 * Purpose:	Tests named data types.
 *
 * Return:	Success:	0
 *
 *		Failure:	number of errors
 *
 * Programmer:	Binh-Minh Ribler (use C version)
 *		January, 2007
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
const H5std_string filename2("dtypes2.h5");
static void test_named ()
{
    static hsize_t	ds_size[2] = {10, 20};
    hsize_t		i;
    unsigned 		attr_data[10][20];
    
    //MESSAGE(5, ("named data types"));
    SUBTEST("Named data types");
    try {
	// Create the file
	H5File file(filename2, H5F_ACC_TRUNC);

	// Create simple dataspace
	DataSpace space(2, ds_size, ds_size);

	// Predefined types cannot be committed
	try {
	    PredType nativeint(PredType::NATIVE_INT);
	    nativeint.commit(file, "test_named_1 (should not exist)");

	    // Should FAIL but didn't, so throw an invalid action exception
	    throw InvalidActionException("PredType::commit", "Attempted to commit a predefined datatype.");
	} catch (DataTypeIException err) {}

	// Copy a predefined data type and commit the copy
        IntType itype(PredType::NATIVE_INT);
        itype.commit(file, "native-int");
	if (itype.committed() <= (bool)0)
	    cerr << "IntType::committed() returned false" << endl;

	// We should not be able to modify a type after it has been committed.
	try {
	    itype.setPrecision(256);

	    // Should FAIL but didn't, so throw an invalid action exception
	    throw InvalidActionException("IntType::setPrecision", "Attempted to modify a committed datatype.");
	} catch (DataTypeIException err) {}

	// We should not be able to re-commit a committed type
	try {
	    itype.commit(file, "test_named_2 (should not exist)");

	    // Should FAIL but didn't, so throw an invalid action exception
	    throw InvalidActionException("IntType::commit", "Attempted to re-commit a committed data type.");
	} catch (DataTypeIException err) {}

	// It should be possible to define an attribute for the named type
	Attribute attr1 = itype.createAttribute("attr1", PredType::NATIVE_UCHAR, space);

	for (i=0; i<ds_size[0]*ds_size[1]; i++)
	    attr_data[0][i] = (int)i; //tricky

	attr1.write(PredType::NATIVE_UINT, attr_data);
	attr1.close();

	// Copying a committed type should result in a transient type which 
	// is not locked.
	IntType trans_type;
	trans_type.copy(itype);
	bool iscommitted = trans_type.committed();
	verify_val(iscommitted, (bool)0, "DataType::committed() - Copying a named type should result in a transient type!", __LINE__, __FILE__);

	trans_type.setPrecision(256);
	trans_type.close();

	// Close the committed type and reopen it.  It should return a 
	// named type.

	/* This had something to do with the way IntType was returned and 
	   assigned and caused itype.committed not working correctly.  So, 
	   use another_type for now.  Need to address it later - BMR 01/2007
	itype = file.openIntType("native-int");
	iscommitted = itype.committed();
	*/
	IntType another_type = file.openIntType("native-int");
	iscommitted = another_type.committed();
	if (!iscommitted)
	    throw InvalidActionException("IntType::committed()", "Opened named types should be named types!");

	// Create a dataset that uses the named type
	DataSet dset = file.createDataSet("dset1", itype, space);

	// Get the dataset's data type and make sure it's a named type
	DataType *ds_type = new DataType(dset.getDataType());
	iscommitted = ds_type->committed();
	if (!iscommitted)
	    throw InvalidActionException("IntType::committed()", "1 Dataset type should be named type!");

	// Close the dataset, then close its type, then reopen the dataset
	dset.close();
	ds_type->close();

	dset = file.openDataSet("dset1");

	// Get the dataset's data type and make sure it's a named type
	ds_type = new DataType(dset.getDataType());
	iscommitted = ds_type->committed();
	if (!iscommitted)
	    throw InvalidActionException("IntType::committed()", "Dataset type should be named type!");

	// Close the dataset and create another with the type returned 
	// from the first dataset.
	dset.close();
	dset = file.createDataSet("dset2", *ds_type, space);

	// Close then reopen the second dataset and make sure the type is 
	// shared
	ds_type->close();
	dset.close();
	dset = file.openDataSet("dset2");
	ds_type = new DataType(dset.getDataType());
	iscommitted = ds_type->committed();
	if (!iscommitted)
	    throw InvalidActionException("IntType::committed()", "Dataset type should be named type!");

	ds_type->close();
    
	// Get the dataset data type by way of copying via the dataset. The
	// result should be modifiable.
	IntType copied_type;
	copied_type.copy(dset);
	copied_type.setPrecision(256);	// modifiable
	copied_type.close();

	// Close all objects and file.
	dset.close();
	itype.close();
	space.close();
	file.close();
    PASSED();
    }   // end of try block
    catch (Exception E) {
        issue_fail_msg(E.getCFuncName(), __LINE__, __FILE__, E.getCDetailMsg());
    }
}   // test_named


/****************************************************************
**
**  test_types(): Main data types testing routine.
**
****************************************************************/
void test_types(void)
{
    // Output message about test being performed
    MESSAGE(5, ("Testing Generic Data Types\n"));

    // Test basic datatypes
    test_classes();
    test_copy();
    test_query();
    test_transient();
    test_named();

}   // test_types()


/*-------------------------------------------------------------------------
 * Function:	cleanup_types
 *
 * Purpose:	Cleanup temporary test files
 *
 * Return:	none
 *
 * Programmer:	Quincey Koziol
 *		September 10, 1999
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
void
cleanup_types(void)
{
    for (int i = 0; i < 3; i++)
	HDremove(FILENAME[i]);
}

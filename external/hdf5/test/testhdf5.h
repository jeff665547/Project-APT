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

/*
 * This header file contains information required for testing the HDF5 library.
 */

#ifndef TESTHDF5_H
#define TESTHDF5_H

/*
 * Include required headers.  This file tests internal library functions,
 * so we include the private headers here.
 */
#include "H5private.h"
#include "H5Eprivate.h"

/* Include generic testing header also */
#include "h5test.h"

/* Use %ld to print the value because long should cover most cases. */
/* Used to make certain a return value _is_not_ a value */
#define CHECK(ret, val, where) do {					      \
    if (GetTestVerbosity()>=VERBO_HI) print_func("   Call to routine: %15s at line %4d " \
				"in %s returned %ld \n",		      \
				where, (int)__LINE__, __FILE__,		      \
				(long)(ret));				      \
    if ((ret) == (val)) {						      \
	TestErrPrintf("*** UNEXPECTED RETURN from %s is %ld at line %4d "     \
		   "in %s\n", where, (long)(ret), (int)__LINE__, __FILE__);   \
	H5Eprint (stdout);						      \
    }									      \
} while(0)

#define CHECK_I(ret,where) {						      \
   if (GetTestVerbosity()>=VERBO_HI) {						      \
      print_func("   Call to routine: %15s at line %4d in %s returned %ld\n", \
                 (where), (int)__LINE__, __FILE__, (long)(ret));	      \
   }									      \
   if ((ret)<0) {							      \
      TestErrPrintf ("*** UNEXPECTED RETURN from %s is %ld line %4d in %s\n", \
                  (where), (long)(ret), (int)__LINE__, __FILE__);	      \
      H5Eprint (stdout);						      \
   }									      \
}

#define CHECK_PTR(ret,where) {						      \
   if (GetTestVerbosity()>=VERBO_HI) {					      \
      print_func("   Call to routine: %15s at line %4d in %s returned %p\n",  \
                 (where), (int)__LINE__, __FILE__, (ret));		      \
   }									      \
   if (!(ret)) {							      \
      TestErrPrintf ("*** UNEXPECTED RETURN from %s is NULL line %4d in %s\n", \
                  (where), (int)__LINE__, __FILE__);			      \
      H5Eprint (stdout);						      \
   }									      \
}

/* Used to make certain a scalar return value _is_ a value */
#define VERIFY(x, val, where) do {					      \
    if (GetTestVerbosity()>=VERBO_HI) {					      \
	print_func("   Call to routine: %15s at line %4d in %s had value "    \
		   "%ld \n", (where), (int)__LINE__, __FILE__, (long)(x));    \
    }									      \
    if ((x) != (val)) {							      \
	TestErrPrintf("*** UNEXPECTED VALUE from %s should be %ld, but is %ld at line %4d " \
		   "in %s\n", (where), (long)(val), (long)(x), (int)__LINE__, __FILE__); \
	H5Eprint (stdout);						      \
    }									      \
} while(0)

/* Used to make certain a string return value _is_ a value */
#define VERIFY_STR(x, val, where) do {					      \
    if (GetTestVerbosity()>=VERBO_HI) {				              \
	print_func("   Call to routine: %15s at line %4d in %s had value "    \
		   "%s \n", (where), (int)__LINE__, __FILE__, x);    \
    }									      \
    if (HDstrcmp(x, val)) {					              \
	TestErrPrintf("*** UNEXPECTED VALUE from %s should be %s, but is %s at line %4d " \
		   "in %s\n", where, val, x, (int)__LINE__, __FILE__);        \
	H5Eprint (stdout);					              \
    }									      \
} while(0)

/* Used to document process through a test and to check for errors */
#define RESULT(ret,func) do {						      \
    if (GetTestVerbosity()>VERBO_MED) {					      \
	print_func("   Call to routine: %15s at line %4d in %s returned "     \
		   "%ld\n", func, (int)__LINE__, __FILE__, (long)(ret));      \
    }									      \
    if (GetTestVerbosity()>=VERBO_HI)					      \
        H5Eprint(stdout);					              \
    if ((ret) == FAIL) {						      \
	TestErrPrintf("*** UNEXPECTED RETURN from %s is %ld at line %4d "     \
		   "in %s\n", func, (long)(ret), (int)__LINE__, __FILE__);    \
	H5Eprint (stdout);						      \
    }									      \
} while(0)

/* Used to document process through a test */
#define MESSAGE(V,A) {if (GetTestVerbosity()>(V)) print_func A;}

/* definitions for command strings */
#define VERBOSITY_STR   "Verbosity"
#define SKIP_STR        "Skip"
#define TEST_STR        "Test"
#define CLEAN_STR       "Cleanup"

#ifdef __cplusplus
extern "C" {
#endif

/* Prototypes for the test routines */
void                    test_metadata(void);
void                    test_tst(void);
void                    test_heap(void);
void                    test_refstr(void);
void                    test_file(void);
void                    test_h5t(void);
void                    test_h5s(void);
void                    test_h5d(void);
void                    test_attr(void);
void                    test_select(void);
void                    test_time(void);
void                    test_reference(void);
void                    test_vltypes(void);
void                    test_vlstrings(void);
void                    test_iterate(void);
void                    test_array(void);
void                    test_genprop(void);
void			test_configure(void);
void			test_misc(void);
void			test_skiplist(void);

/* Prototypes for the cleanup routines */
void                    cleanup_metadata(void);
void                    cleanup_file(void);
void                    cleanup_h5s(void);
void                    cleanup_attr(void);
void                    cleanup_select(void);
void                    cleanup_time(void);
void                    cleanup_reference(void);
void                    cleanup_vltypes(void);
void                    cleanup_vlstrings(void);
void                    cleanup_iterate(void);
void                    cleanup_array(void);
void                    cleanup_genprop(void);
void			cleanup_configure(void);
void			cleanup_misc(void);

#ifdef __cplusplus
}
#endif
#endif /* TESTHDF5_H */

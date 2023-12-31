#! /bin/sh
#
# Copyright by The HDF Group.
# Copyright by the Board of Trustees of the University of Illinois.
# All rights reserved.
#
# This file is part of HDF5.  The full HDF5 copyright notice, including
# terms governing use, modification, and redistribution, is contained in
# the files COPYING and Copyright.html.  COPYING can be found at the root
# of the source code distribution tree; Copyright.html can be found at the
# root level of an installed copy of the electronic HDF5 document set and
# is linked from the top-level documents page.  It can also be found at
# http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have
# access to either file, you may request a copy from help@hdfgroup.org.
#
# Tests for the h5ls tool

H5LS=h5ls               # The tool name
H5LS_BIN=`pwd`/$H5LS    # The path of the tool binary

CMP='cmp -s'
DIFF='diff -c'
NLINES=20			# Max. lines of output to display if test fails

nerrors=0
verbose=yes

# The build (current) directory might be different than the source directory.
if test -z "$srcdir"; then
    srcdir=.
fi
test -d ../testfiles || mkdir ../testfiles

# Print a line-line message left justified in a field of 70 characters
# beginning with the word "Testing".
TESTING() {
    SPACES="                                                               "
    echo "Testing $* $SPACES" |cut -c1-70 |tr -d '\012'
}

# Run a test and print PASS or *FAIL*. For now, if h5ls can complete
# with exit status 0, consider it pass. If a test fails then increment
# the `nerrors' global variable and (if $verbose is set) display up to $NLINS
# lines of the actual output from the tool test.  The actual output is not
# removed if $HDF5_NOCLEANUP has a non-zero value.
# Arguemnts:
# $1 -- actual output filename to use
# $2 and on -- argument for the h5ls tool
TOOLTEST() {
    expect="$srcdir/../testfiles/$1"
    actual="../testfiles/`basename $1 .ls`.out"
    shift

    # Run test.
    # Stderr is included in stdout so that the diff can detect
    # any unexpected output from that stream too.
    TESTING $H5LS $@
    (
	echo "#############################"
	echo " output for '$H5LS $@'" 
	echo "#############################"
	cd $srcdir/../testfiles
        $RUNSERIAL $H5LS_BIN "$@"
    ) 2>&1 |sed 's/Modified:.*/Modified:  XXXX-XX-XX XX:XX:XX XXX/' >$actual
    
    exitcode=$?
    if [ $exitcode -ne 0 ]; then
	echo "*FAILED*"
	nerrors="`expr $nerrors + 1`"
	if [ yes = "$verbose" ]; then
	    echo "test returned with exit code $exitcode"
	    echo "test output: (up to $NLINES lines)"
	    head -$NLINES $actual
	    echo "***end of test output***"
	    echo ""
	fi
    elif [ ! -f $expect ]; then
	# Create the expect file if it doesn't yet exist.
        echo " CREATED"
        cp $actual $expect
    elif $CMP $expect $actual; then
        echo " PASSED"
    else
        echo "*FAILED*"
	echo "    Expected result differs from actual result"
	nerrors="`expr $nerrors + 1`"
	test yes = "$verbose" && $DIFF $expect $actual |sed 's/^/    /'
    fi

    # Clean up output file
    if test -z "$HDF5_NOCLEANUP"; then
	rm -f $actual
    fi
}

##############################################################################
##############################################################################
###			  T H E   T E S T S                                ###
##############################################################################
##############################################################################

# Toss in a bunch of tests.  Not sure if they are the right kinds.
# test the help syntax
TOOLTEST help-1.ls -w80 -h
TOOLTEST help-2.ls -w80 -help
TOOLTEST help-3.ls -w80 -?

# test simple command
TOOLTEST tall-1.ls -w80 tall.h5
TOOLTEST tall-2.ls -w80 -r -d tall.h5
TOOLTEST tgroup.ls -w80 tgroup.h5

# test for displaying groups
TOOLTEST tgroup-1.ls -w80 -r -g tgroup.h5

# test for displaying simple space datasets
TOOLTEST tdset-1.ls -w80 -r -d tdset.h5

# test for displaying soft links
TOOLTEST tslink-1.ls -w80 -r tslink.h5

# tests for hard links
TOOLTEST thlink-1.ls -w80 thlink.h5

# tests for compound data types
TOOLTEST tcomp-1.ls -w80 -r -d tcompound.h5

#test for the nested compound type
TOOLTEST tnestcomp-1.ls -w80 -r -d tnestedcomp.h5

# test for loop detection
TOOLTEST tloop-1.ls -w80 -r -d tloop.h5

# test for string 
TOOLTEST tstr-1.ls -w80 -r -d tstr.h5

# test test file created from lib SAF team
TOOLTEST tsaf.ls -w80 -r -d tsaf.h5

# test for variable length data types
TOOLTEST tvldtypes1.ls -w80 -r -d tvldtypes1.h5

# test for array data types
TOOLTEST tarray1.ls -w80 -r -d tarray1.h5

# test for empty data
TOOLTEST tempty.ls -w80 -d tempty.h5

# test for all dataset types written to attributes
# enable -S for avoiding printing NATIVE types
TOOLTEST tattr2.ls -w80 -v -S tattr2.h5

# tests for error handling.
# test for non-existing file
TOOLTEST nosuchfile.ls nosuchfile.h5

if test $nerrors -eq 0 ; then
	echo "All h5ls tests passed."
fi

exit $nerrors

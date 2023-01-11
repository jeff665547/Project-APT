#
# # cvs:affy/sdk/copynumber/apt-copynumber-familial/regression/Makefile ---
#
# $Id: Makefile,v 1.3 2009-11-02 22:56:14 jbiela Exp $
#

t_pgms:=test-regression-copynumber-familial test-regression-copynumber-familial-quick
a_pgms:=${t_pgms} apt-copynumber-familial

#
sdk_root:=../../..
include ${sdk_root}/Makefile.defs
#
sdk_exe+=${t_pgms}
_regression_run+=cp_pgms regression_1
_integration_run+=cp_pgms integration_1
# _valgrind_run+=cp_pgms valgrind_1
_clean_run+=clean_regression
#
include ${sdk_makefile_post}

cp_pgms:
	${CP} ${a_pgms:%=${sdk_output_bin}/%} .

#valgrind_1:
#	./test-regression-copynumber-familial-quick.valgrind
#
integration_1:
	./test-regression-copynumber-familial-quick

regression_1:
	./test-regression-copynumber-familial

clean_regression:
	-rm -r test-generated


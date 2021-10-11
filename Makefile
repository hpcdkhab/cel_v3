#     
# File:  Makefile
# Project CRESTA (see details on https://cresta-project.eu) Exascale library
# Send you email about the bugs to
# The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
# Dmitry Khabi,  (email: khabi@hlrs.de)
# Created on March 25, 2013
#

export PROJECT_DIR=$(shell pwd)

#export TASK4_LIB_DIR=$(PROJECT_DIR)/ext_lib/lib_hlrs_crayxe6_cray_debug_double_i64
#export TASK4_MOD_DIR=$(PROJECT_DIR)/ext_lib/mod/mod_hlrs_crayxe6_cray_debug_double_i64

export TASK4_LIB_DIR=$(PROJECT_DIR)/ext_lib/
export TASK4_MOD_DIR=$(PROJECT_DIR)/ext_lib/mod/

export ROOT_OBJ_DIR = $(PROJECT_DIR)/obj
export ROOT_INC_DIR = $(PROJECT_DIR)/src/
export ROOT_MOD_DIR= $(PROJECT_DIR)/src/
export LIB_DIR = $(PROJECT_DIR)/lib

include ./make_env.inc
include ./makefile_base.inc
include ./makefile_fbase.inc
include ./makefile_fcomm.inc
include ./makefile_fspmt.inc
include ./makefile_fbasedistr.inc
include ./makefile_fspmtconverter.inc
include ./makefile_fspmtdistr.inc
include ./makefile_fprec.inc
include ./makefile_fop.inc
include ./makefile_fcgalg.inc
include ./makefile_felmsolver.inc


.PHONY: all tests

all: base fbase fcomm fspmt fbasedistr fspmtconverter fspmtdistr fprec fop fcgalg felmsolver
default:  base fbase fcomm fspmt fbasedistr fspmtconverter fspmtdistr fprec fop fcgalg felmsolver

fort: fbase fcomm fspmt fbasedistr fspmtconverter fspmtdistr fprec fop fcgalg felmsolver

$(LIB_DIR) :
	mkdir $(LIB_DIR) -p




tests:
	$(MAKE) -C tests all

clean_tests:
	$(MAKE) -C tests clean


clean_list=clean_base clean_fbase clean_fcomm clean_fspmt clean_fbasedistr clean_fspmtconverter clean_fspmtdistr clean_fprec clean_fop clean_fcgalg clean_felmsolver

clean: $(clean_list) 
	rm -rf $(ROOT_OBJ_DIR) $(LIB_DIR)
	rm -rf $(COMPILER_DIR)/cray.pl
	$(MAKE) -C tests clean

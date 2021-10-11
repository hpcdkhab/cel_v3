#     
# File:  Makefile for library libcelfprec.a
# Project CRESTA Exascale library (see details on https://www.cresta-project.eu)
# Send you email about the bugs to
# The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
# Dmitry Khabi,  (email: khabi@hlrs.de)
# Created on March 25, 2013
#


# source files.
FPREC_SRC = cel_prec_module.f90
#output library 
export FPREC_OUT = $(LIB_DIR)/libcelfprec.a

#directories

export FPREC_SRC_DIR = $(PROJECT_DIR)/src/fprec
export FPREC_OBJ_DIR = $(ROOT_OBJ_DIR)/fprec
export FPREC_MOD_DIR = $(ROOT_MOD_DIR)/include/fprec
export FPREC_INC_DIR = $(ROOT_INC_DIR)/include/fprec

.PHONY: FPREC_directories
# include directories
FPREC_INCLUDES_TMP = $(FPREC_INC_DIR) $(FSPMTDISTR_INC_DIR)  $(FSPMTCONVERTER_INC_DIR) $(FCOMM_INC_DIR) $(FSPMT_INC_DIR) $(FBASE_INC_DIR)  
export FPREC_INCLUDES=$(foreach d, $(FPREC_INCLUDES_TMP), -I$d)

FPREC_SOURCES=$(foreach file,$(FPREC_SRC),$(FPREC_SRC_DIR)/$(file))
FPREC_OBJECTS=$(foreach file,$(FPREC_SRC),$(FPREC_OBJ_DIR)/$(basename $(file)).o)


fprec: FPREC_directories $(FPREC_OUT)

$(FPREC_OBJ_DIR)/%.o:	$(FPREC_SRC_DIR)/%.f90 
	$(MPI_F90) $(FPREC_INCLUDES) $(FFLAGS) $(F90_MODULE_DIR) $(FPREC_MOD_DIR) -c $< -o $@

$(FPREC_OUT): $(FPREC_OBJECTS)  $(LIB_DIR) 
	ar rcs $(FPREC_OUT) $(FPREC_OBJECTS)

FPREC_directories: $(FPREC_MOD_DIR) $(FPREC_OBJ_DIR)

$(FPREC_OBJ_DIR): 
	mkdir $(FPREC_OBJ_DIR) -p
$(FPREC_MOD_DIR): 
	mkdir $(FPREC_MOD_DIR) -p


clean_fprec:
	rm -rf $(FPREC_OBJECTS) $(FPREC_OUT)
	rm -rf $(FPREC_MOD_DIR)


#     
# File:  Makefile for FCOMM directory, library libcelfcomm.a
# Project CRESTA (see details on https://www.cresta-project.eu) Exascale library
# Send you email about the bugs to
# The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
# Dmitry Khabi,  (email: khabi@hlrs.de)
# Created on Nov 13, 2013
#


# source files.
FCOMM_SRC = cel_comm_module.f90


#output library 
export FCOMM_OUT = $(LIB_DIR)/libcelfcomm.a

#directories

export FCOMM_SRC_DIR = $(PROJECT_DIR)/src/fcomm
export FCOMM_OBJ_DIR = $(ROOT_OBJ_DIR)/fcomm
export FCOMM_MOD_DIR = $(ROOT_MOD_DIR)/include/fcomm
export FCOMM_INC_DIR = $(ROOT_INC_DIR)/include/fcomm

.PHONY: FCOMM_directories

#include directories
FCOMM_INCLUDES_TMP = $(FCOMM_INC_DIR) $(FBASE_INC_DIR) 
export FCOMM_INCLUDES=$(foreach d, $(FCOMM_INCLUDES_TMP), -I$d)

FCOMM_SOURCES=$(foreach file,$(FCOMM_SRC),$(FCOMM_SRC_DIR)/$(file))
FCOMM_OBJECTS=$(foreach file,$(FCOMM_SRC),$(FCOMM_OBJ_DIR)/$(basename $(file)).o)


fcomm: FCOMM_directories $(FCOMM_OUT)

$(FCOMM_OBJ_DIR)/%.o:	$(FCOMM_SRC_DIR)/%.f90 
	$(MPI_F90) $(FCOMM_INCLUDES) $(FFLAGS) $(F90_MODULE_DIR) $(FCOMM_MOD_DIR) -c $< -o $@

$(FCOMM_OUT): $(FCOMM_OBJECTS)  $(LIB_DIR) 
	ar rcs $(FCOMM_OUT) $(FCOMM_OBJECTS)

FCOMM_directories: $(FCOMM_MOD_DIR) $(FCOMM_OBJ_DIR)

$(FCOMM_OBJ_DIR): 
	mkdir $(FCOMM_OBJ_DIR) -p
$(FCOMM_MOD_DIR): 
	mkdir $(FCOMM_MOD_DIR) -p


clean_fcomm:
	rm -rf $(FCOMM_OBJECTS) $(FCOMM_OUT)
	rm -rf $(FCOMM_MOD_DIR)


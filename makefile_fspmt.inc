#     
# File:  Makefile for FSPMT directory, library libcelFSPMT.a
# Project CRESTA (see details on https://cresta-project.eu) Exascale library
# Send you email about the bugs to
# The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
# Dmitry Khabi,  (email: khabi@hlrs.de)
# Created on March 25, 2013
#


# source files.
FSPMT_SRC = cel_vec_module.f90
#FSPMT_SRC += cel_cg_module.f90
FSPMT_SRC += cel_domain_3D_module.f90
FSPMT_SRC += cel_discretisation_task_module.f90
FSPMT_SRC += cel_sp_mat_module.f90
FSPMT_SRC += cel_c_sp_mat_interface_module.f90
FSPMT_SRC += cel_sp_mat_generate_module.f90
FSPMT_SRC += cel_sp_mat_check_module.f90

#output library 
export FSPMT_OUT = $(LIB_DIR)/libcelfspmt.a

#directories

export FSPMT_SRC_DIR = $(PROJECT_DIR)/src/fspmt
export FSPMT_OBJ_DIR = $(ROOT_OBJ_DIR)/fspmt
export FSPMT_MOD_DIR = $(ROOT_MOD_DIR)/include/fspmt
export FSPMT_INC_DIR = $(ROOT_INC_DIR)/include/fspmt

.PHONY: fspmt_directories

# include directories
FSPMT_INCLUDES_TMP = $(FSPMT_INC_DIR) $(FBASE_INC_DIR) 
export FSPMT_INCLUDES=$(foreach d, $(FSPMT_INCLUDES_TMP), -I$d)

FSPMT_SOURCES=$(foreach file,$(FSPMT_SRC),$(FSPMT_SRC_DIR)/$(file))
FSPMT_OBJECTS=$(foreach file,$(FSPMT_SRC),$(FSPMT_OBJ_DIR)/$(basename $(file)).o)


fspmt: fspmt_directories $(FSPMT_OUT)

$(FSPMT_OBJ_DIR)/%.o:	$(FSPMT_SRC_DIR)/%.f90 
	$(F90) $(FSPMT_INCLUDES) $(FFLAGS) $(F90_MODULE_DIR) $(FSPMT_MOD_DIR) -c $< -o $@

$(FSPMT_OUT): $(FSPMT_OBJECTS)  $(LIB_DIR) 
	ar rcs $(FSPMT_OUT) $(FSPMT_OBJECTS)

fspmt_directories: $(FSPMT_MOD_DIR) $(FSPMT_OBJ_DIR)

$(FSPMT_OBJ_DIR): 
	mkdir $(FSPMT_OBJ_DIR) -p
$(FSPMT_MOD_DIR): 
	mkdir $(FSPMT_MOD_DIR) -p


clean_fspmt:
	rm -rf $(FSPMT_OBJECTS) $(FSPMT_OUT)
	rm -rf $(FSPMT_MOD_DIR)

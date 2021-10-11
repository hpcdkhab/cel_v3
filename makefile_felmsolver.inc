#     
# File:  Makefile for FSPMT directory, library libcelFSPMT.a
# Project CRESTA (see details on https://www.cresta-project.eu) Exascale library
# Send you email about the bugs to
# The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
# Dmitry Khabi,  (email: khabi@hlrs.de)
# Created on March 25, 2013
#


# source files.
FELMSOLVER_SRC = cel_elmsolver_module.f90


#output library 
export FELMSOLVER_OUT = $(LIB_DIR)/libcelfelmsolver.a

#directories

export FELMSOLVER_SRC_DIR = $(PROJECT_DIR)/src/felmsolver
export FELMSOLVER_OBJ_DIR = $(ROOT_OBJ_DIR)/felmsolver
export FELMSOLVER_MOD_DIR = $(ROOT_MOD_DIR)/include/felmsolver
export FELMSOLVER_INC_DIR = $(ROOT_INC_DIR)/include/felmsolver

.PHONY: FELMSOLVER_directories

# include directories
FELMSOLVER_INCLUDES_TMP = $(FELMSOLVER_INC_DIR) $(FPREC_INC_DIR) $(FSPMTDISTR_INC_DIR) $(FSPMTCONVERTER_INC_DIR) $(FCOMM_INC_DIR) $(FSPMT_INC_DIR) $(FBASEDISTR_INC_DIR) $(FBASE_INC_DIR) $(FOP_INC_DIR) $(FCGALG_INC_DIR)
export FELMSOLVER_INCLUDES=$(foreach d, $(FELMSOLVER_INCLUDES_TMP), -I$d)

FELMSOLVER_SOURCES=$(foreach file,$(FELMSOLVER_SRC),$(FELMSOLVER_SRC_DIR)/$(file))
FELMSOLVER_OBJECTS=$(foreach file,$(FELMSOLVER_SRC),$(FELMSOLVER_OBJ_DIR)/$(basename $(file)).o)


felmsolver: felmsolver_directories $(FELMSOLVER_OUT)

$(FELMSOLVER_OBJ_DIR)/%.o:	$(FELMSOLVER_SRC_DIR)/%.f90 
	$(MPI_F90) $(FELMSOLVER_INCLUDES) $(FFLAGS) $(F90_MODULE_DIR) $(FELMSOLVER_MOD_DIR) -c $< -o $@

$(FELMSOLVER_OUT): $(FELMSOLVER_OBJECTS)  $(LIB_DIR) 
	ar rcs $(FELMSOLVER_OUT) $(FELMSOLVER_OBJECTS)

felmsolver_directories: $(FELMSOLVER_MOD_DIR) $(FELMSOLVER_OBJ_DIR)

$(FELMSOLVER_OBJ_DIR): 
	mkdir $(FELMSOLVER_OBJ_DIR) -p
$(FELMSOLVER_MOD_DIR): 
	mkdir $(FELMSOLVER_MOD_DIR) -p


clean_felmsolver:
	rm -rf $(FELMSOLVER_OBJECTS) $(FELMSOLVER_OUT)
	rm -rf $(FELMSOLVER_MOD_DIR)


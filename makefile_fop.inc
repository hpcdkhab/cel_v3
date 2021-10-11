#     
# File:  Makefile for FSPMT directory, library libcelFSPMT.a
# Project CRESTA (see details on https://cresta-project.eu) Exascale library
# Send you email about the bugs to
# The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
# Dmitry Khabi,  (email: khabi@hlrs.de)
# Created on March 25, 2013
#


# source files.
FOP_SRC = cel_fop_mv_module.f90
FOP_SRC += cel_fop_vec_distr_module.f90
FOP_SRC += cel_fop_block_jacobi_module.f90

#output library 
export FOP_OUT = $(LIB_DIR)/libcelfop.a

#directories

export FOP_SRC_DIR = $(PROJECT_DIR)/src/fop
export FOP_OBJ_DIR = $(ROOT_OBJ_DIR)/fop
export FOP_MOD_DIR = $(ROOT_MOD_DIR)/include/fop
export FOP_INC_DIR = $(ROOT_INC_DIR)/include/fop

.PHONY: FOP_directories

# include directories
FOP_INCLUDES_TMP = $(FOP_INC_DIR) $(FSPMTDISTR_INC_DIR) $(FCOMM_INC_DIR) $(FSPMT_INC_DIR) $(FBASE_INC_DIR) 
export FOP_INCLUDES=$(foreach d, $(FOP_INCLUDES_TMP), -I$d)

FOP_SOURCES=$(foreach file,$(FOP_SRC),$(FOP_SRC_DIR)/$(file))
FOP_OBJECTS=$(foreach file,$(FOP_SRC),$(FOP_OBJ_DIR)/$(basename $(file)).o)


fop: FOP_directories $(FOP_OUT)

$(FOP_OBJ_DIR)/%.o:	$(FOP_SRC_DIR)/%.f90 
	$(MPI_F90) $(FOP_INCLUDES) $(FFLAGS) $(F90_MODULE_DIR) $(FOP_MOD_DIR) -c $< -o $@

$(FOP_OUT): $(FOP_OBJECTS)  $(LIB_DIR) 
	ar rcs $(FOP_OUT) $(FOP_OBJECTS)

FOP_directories: $(FOP_MOD_DIR) $(FOP_OBJ_DIR)

$(FOP_OBJ_DIR): 
	mkdir $(FOP_OBJ_DIR) -p
$(FOP_MOD_DIR): 
	mkdir $(FOP_MOD_DIR) -p


clean_fop:
	rm -rf $(FOP_OBJECTS) $(FOP_OUT)
	rm -rf $(FOP_MOD_DIR)


#     
# File:  Makefile for FSPMT directory, library libcelFSPMT.a
# Project CRESTA (see details on https://www.cresta-project.eu) Exascale library
# Send you email about the bugs to
# The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
# Dmitry Khabi,  (email: khabi@hlrs.de)
# Created on March 25, 2013
#


# source files.
FCGALG_SRC = cel_alg_module.f90
FCGALG_SRC += cel_cgalg_module.f90
FCGALG_SRC += cel_gmresalg_module.f90

#output library 
export FCGALG_OUT = $(LIB_DIR)/libcelcgalg.a

#directories

export FCGALG_SRC_DIR = $(PROJECT_DIR)/src/fcgalg
export FCGALG_OBJ_DIR = $(ROOT_OBJ_DIR)/fcgalg
export FCGALG_MOD_DIR = $(ROOT_MOD_DIR)/include/fcgalg
export FCGALG_INC_DIR = $(ROOT_INC_DIR)/include/fcgalg

.PHONY: FCGALG_directories

# include directories
FCGALG_INCLUDES_TMP = $(FCGALG_INC_DIR) $(FSPMTDISTR_INC_DIR) $(FCOMM_INC_DIR) $(FSPMT_INC_DIR) $(FBASE_INC_DIR) $(FOP_INC_DIR)
export FCGALG_INCLUDES=$(foreach d, $(FCGALG_INCLUDES_TMP), -I$d)

FCGALG_SOURCES=$(foreach file,$(FCGALG_SRC),$(FCGALG_SRC_DIR)/$(file))
FCGALG_OBJECTS=$(foreach file,$(FCGALG_SRC),$(FCGALG_OBJ_DIR)/$(basename $(file)).o)


fcgalg: FCGALG_directories $(FCGALG_OUT)

$(FCGALG_OBJ_DIR)/%.o:	$(FCGALG_SRC_DIR)/%.f90 
	$(MPI_F90) $(FCGALG_INCLUDES) $(FFLAGS) $(F90_MODULE_DIR) $(FCGALG_MOD_DIR) -c $< -o $@

$(FCGALG_OUT): $(FCGALG_OBJECTS)  $(LIB_DIR) 
	ar rcs $(FCGALG_OUT) $(FCGALG_OBJECTS)

FCGALG_directories: $(FCGALG_MOD_DIR) $(FCGALG_OBJ_DIR)

$(FCGALG_OBJ_DIR): 
	mkdir $(FCGALG_OBJ_DIR) -p
$(FCGALG_MOD_DIR): 
	mkdir $(FCGALG_MOD_DIR) -p


clean_fcgalg:
	rm -rf $(FCGALG_OBJECTS) $(FCGALG_OUT)
	rm -rf $(FCGALG_MOD_DIR)


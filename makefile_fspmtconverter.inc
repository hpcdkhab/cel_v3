#     
# File:  Makefile for library libcelfspmtconverter.a
# Project CRESTA Exascale library (see details on https://www.cresta-project.eu)
# Send you email about the bugs to
# The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
# Dmitry Khabi,  (email: khabi@hlrs.de)
# Created on March 25, 2013
#


# source files.
FSPMTCONVERTER_SRC = cel_sp_mat_jad_converter_module.f90



#output library 
export FSPMTCONVERTER_OUT = $(LIB_DIR)/libcelfspmtconverter.a

#directories

export FSPMTCONVERTER_SRC_DIR = $(PROJECT_DIR)/src/fspmtconverter
export FSPMTCONVERTER_OBJ_DIR = $(ROOT_OBJ_DIR)/fspmtconverter
export FSPMTCONVERTER_MOD_DIR = $(ROOT_MOD_DIR)/include/fspmtconverter
export FSPMTCONVERTER_INC_DIR = $(ROOT_INC_DIR)/include/fspmtconverter

.PHONY: FSPMTCONVERTER_directories

# include directories
FSPMTCONVERTER_INCLUDES_TMP = $(FSPMTCONVERTER_INC_DIR) $(FCOMM_INC_DIR) $(FSPMT_INC_DIR) $(FBASE_INC_DIR) 
export FSPMTCONVERTER_INCLUDES=$(foreach d, $(FSPMTCONVERTER_INCLUDES_TMP), -I$d)

FSPMTCONVERTER_SOURCES=$(foreach file,$(FSPMTCONVERTER_SRC),$(FSPMTCONVERTER_SRC_DIR)/$(file))
FSPMTCONVERTER_OBJECTS=$(foreach file,$(FSPMTCONVERTER_SRC),$(FSPMTCONVERTER_OBJ_DIR)/$(basename $(file)).o)


fspmtconverter: FSPMTCONVERTER_directories $(FSPMTCONVERTER_OUT)

$(FSPMTCONVERTER_OBJ_DIR)/%.o:	$(FSPMTCONVERTER_SRC_DIR)/%.f90 
	$(MPI_F90) $(FSPMTCONVERTER_INCLUDES) $(FFLAGS) $(F90_MODULE_DIR) $(FSPMTCONVERTER_MOD_DIR) -c $< -o $@

$(FSPMTCONVERTER_OUT): $(FSPMTCONVERTER_OBJECTS)  $(LIB_DIR) 
	ar rcs $(FSPMTCONVERTER_OUT) $(FSPMTCONVERTER_OBJECTS)

FSPMTCONVERTER_directories: $(FSPMTCONVERTER_MOD_DIR) $(FSPMTCONVERTER_OBJ_DIR)

$(FSPMTCONVERTER_OBJ_DIR): 
	mkdir $(FSPMTCONVERTER_OBJ_DIR) -p
$(FSPMTCONVERTER_MOD_DIR): 
	mkdir $(FSPMTCONVERTER_MOD_DIR) -p


clean_fspmtconverter:
	rm -rf $(FSPMTCONVERTER_OBJECTS) $(FSPMTCONVERTER_OUT)
	rm -rf $(FSPMTCONVERTER_MOD_DIR)


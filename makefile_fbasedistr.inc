#     
# File:  Makefile for library libcelfspmtdistr.a
# Project CRESTA Exascale library (see details on https://www.cresta-project.eu)
# Send you email about the bugs to
# The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
# Dmitry Khabi,  (email: khabi@hlrs.de)
# Created on March 25, 2013
#


# source files.
FBASEDISTR_SRC = cel_perf_distr_module.f90



#output library 
export FBASEDISTR_OUT = $(LIB_DIR)/libcelfbasedistr.a

#directories

export FBASEDISTR_SRC_DIR = $(PROJECT_DIR)/src/fbasedistr
export FBASEDISTR_OBJ_DIR = $(ROOT_OBJ_DIR)/fbasedistr
export FBASEDISTR_MOD_DIR = $(ROOT_MOD_DIR)/include/fbasedistr
export FBASEDISTR_INC_DIR = $(ROOT_INC_DIR)/include/fbasedistr

.PHONY: FBASEDISTR_directories

# include directories
FBASEDISTR_INCLUDES_TMP = $(FBASEDISTR_INC_DIR) $(FCOMM_INC_DIR) $(FBASE_INC_DIR) 
export FBASEDISTR_INCLUDES=$(foreach d, $(FBASEDISTR_INCLUDES_TMP), -I$d)

FBASEDISTR_SOURCES=$(foreach file,$(FBASEDISTR_SRC),$(FBASEDISTR_SRC_DIR)/$(file))
FBASEDISTR_OBJECTS=$(foreach file,$(FBASEDISTR_SRC),$(FBASEDISTR_OBJ_DIR)/$(basename $(file)).o)


fbasedistr: FBASEDISTR_directories $(FBASEDISTR_OUT)

$(FBASEDISTR_OBJ_DIR)/%.o:	$(FBASEDISTR_SRC_DIR)/%.f90 
	$(MPI_F90) $(FBASEDISTR_INCLUDES) $(FFLAGS) $(F90_MODULE_DIR) $(FBASEDISTR_MOD_DIR) -c $< -o $@

$(FBASEDISTR_OUT): $(FBASEDISTR_OBJECTS)  $(LIB_DIR) 
	ar rcs $(FBASEDISTR_OUT) $(FBASEDISTR_OBJECTS)

FBASEDISTR_directories: $(FBASEDISTR_MOD_DIR) $(FBASEDISTR_OBJ_DIR)

$(FBASEDISTR_OBJ_DIR): 
	mkdir $(FBASEDISTR_OBJ_DIR) -p
$(FBASEDISTR_MOD_DIR): 
	mkdir $(FBASEDISTR_MOD_DIR) -p


clean_fbasedistr:
	rm -rf $(FBASEDISTR_OBJECTS) $(FBASEDISTR_OUT)
	rm -rf $(FBASEDISTR_MOD_DIR)


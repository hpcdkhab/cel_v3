#     
# File:  Makefile for library libcelfspmtdistr.a
# Project CRESTA Exascale library (see details on https://www.cresta-project.eu)
# Send you email about the bugs to
# The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
# Dmitry Khabi,  (email: khabi@hlrs.de)
# Created on March 25, 2013
#


# source files.
FSPMTDISTR_SRC = cel_sp_mat_distr_opt_module.f90
FSPMTDISTR_SRC += cel_sp_mat_distr_gl_module.f90
FSPMTDISTR_SRC += cel_sp_mat_distr_module.f90
FSPMTDISTR_SRC += cel_sp_mat_distr_vec_module.f90
FSPMTDISTR_SRC += cel_sp_mat_distr_format_module.f90
FSPMTDISTR_SRC += cel_sp_mat_distr_check_module.f90
#output library 
export FSPMTDISTR_OUT = $(LIB_DIR)/libcelfspmtdistr.a

#directories

export FSPMTDISTR_SRC_DIR = $(PROJECT_DIR)/src/fspmtdistr
export FSPMTDISTR_OBJ_DIR = $(ROOT_OBJ_DIR)/fspmtdistr
export FSPMTDISTR_MOD_DIR = $(ROOT_MOD_DIR)/include/fspmtdistr
export FSPMTDISTR_INC_DIR = $(ROOT_INC_DIR)/include/fspmtdistr

.PHONY: FSPMTDISTR_directories
# include directories
FSPMTDISTR_INCLUDES_TMP = $(FSPMTDISTR_INC_DIR) $(FSPMTCONVERTER_INC_DIR) $(FCOMM_INC_DIR) $(FSPMT_INC_DIR) $(FBASE_INC_DIR) 
export FSPMTDISTR_INCLUDES=$(foreach d, $(FSPMTDISTR_INCLUDES_TMP), -I$d)

FSPMTDISTR_SOURCES=$(foreach file,$(FSPMTDISTR_SRC),$(FSPMTDISTR_SRC_DIR)/$(file))
FSPMTDISTR_OBJECTS=$(foreach file,$(FSPMTDISTR_SRC),$(FSPMTDISTR_OBJ_DIR)/$(basename $(file)).o)


fspmtdistr: FSPMTDISTR_directories $(FSPMTDISTR_OUT)

$(FSPMTDISTR_OBJ_DIR)/%.o:	$(FSPMTDISTR_SRC_DIR)/%.f90 
	$(MPI_F90) $(FSPMTDISTR_INCLUDES) $(FFLAGS) $(F90_MODULE_DIR) $(FSPMTDISTR_MOD_DIR) -c $< -o $@

$(FSPMTDISTR_OUT): $(FSPMTDISTR_OBJECTS)  $(LIB_DIR) 
	ar rcs $(FSPMTDISTR_OUT) $(FSPMTDISTR_OBJECTS)

FSPMTDISTR_directories: $(FSPMTDISTR_MOD_DIR) $(FSPMTDISTR_OBJ_DIR)

$(FSPMTDISTR_OBJ_DIR): 
	mkdir $(FSPMTDISTR_OBJ_DIR) -p
$(FSPMTDISTR_MOD_DIR): 
	mkdir $(FSPMTDISTR_MOD_DIR) -p


clean_fspmtdistr:
	rm -rf $(FSPMTDISTR_OBJECTS) $(FSPMTDISTR_OUT)
	rm -rf $(FSPMTDISTR_MOD_DIR)


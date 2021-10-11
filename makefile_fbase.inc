#     
# File:  Makefile for fbase directory, library libcelfbase.a
# Project CRESTA (see details on https://cresta-project.eu) Exascale library
# Send you email about the bugs to
# The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
# Dmitry Khabi,  (email: khabi@hlrs.de)
# Created on March 25, 2013
#


# source files.
FBASE_SRC = cel_types_module.f90 
FBASE_SRC += cel_cpu_param_module.f90 
FBASE_SRC += cel_mmio_module
FBASE_SRC += cel_algorithms_param_module.f90
FBASE_SRC += cel_timer_module
FBASE_SRC += cel_timer_interface_module
FBASE_SRC += cel_error_module.f90
FBASE_SRC += cel_string_module.f90
FBASE_SRC += cel_interaction_module.f90
FBASE_SRC += cel_omp_module.f90
FBASE_SRC += cel_perf_module.f90
FBASE_SRC += cel_profile_module.f90
FBASE_SRC += cel_extract_command_module.f90
FBASE_SRC += cel_base_print_module.f90
FBASE_SRC += cel_input_parameters_module.f90
FBASE_SRC += cel_vtkio_module.f90
FBASE_SRC += cel_omp_shared_work_module.f90






#output library 
export FBASE_OUT = $(LIB_DIR)/libcelfbase.a

#directories

export FBASE_SRC_DIR = $(PROJECT_DIR)/src/fbase
export FBASE_OBJ_DIR = $(ROOT_OBJ_DIR)/fbase
export FBASE_MOD_DIR = $(ROOT_MOD_DIR)/include/fbase
export FBASE_INC_DIR = $(ROOT_INC_DIR)/include/fbase

.PHONY: fbase_directories

# include directories
FBASE_INCLUDES_TMP = $(FBASE_INC_DIR)  
export FBASE_INCLUDES=$(foreach d, $(FBASE_INCLUDES_TMP), -I$d)

FBASE_SOURCES=$(foreach file,$(FBASE_SRC),$(FBASE_SRC_DIR)/$(file))
FBASE_OBJECTS=$(foreach file,$(FBASE_SRC),$(FBASE_OBJ_DIR)/$(basename $(file)).o)


fbase: $(FBASE_OUT)

$(FBASE_OBJ_DIR)/%.o:	$(FBASE_SRC_DIR)/%.f90 $(FBASE_OBJ_DIR)
	$(F90) $(FBASE_INCLUDES) $(FFLAGS)  $(F90_MODULE_DIR) $(FBASE_MOD_DIR) -c $< -o $@

$(FBASE_OUT): fbase_directories $(FBASE_OBJECTS) 
	ar rcs $(FBASE_OUT) $(FBASE_OBJECTS)

fbase_directories: $(FBASE_OBJ_DIR) $(FBASE_MOD_DIR)  $(LIB_DIR) 

$(FBASE_OBJ_DIR) :
	mkdir $(FBASE_OBJ_DIR) -p
$(FBASE_MOD_DIR) :
	mkdir $(FBASE_MOD_DIR) -p

clean_fbase:
	rm -rf $(FBASE_OBJECTS) $(FBASE_OUT)
	rm -rf $(FBASE_MOD_DIR)

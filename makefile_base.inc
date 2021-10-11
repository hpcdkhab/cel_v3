#     
# File:  Makefile for base directory, library libcelbase.a
# Project CRESTA (see details on https://cresta-project.eu) Exascale library
# Send you email about the bugs to
# The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
# Dmitry Khabi,  (email: khabi@hlrs.de)
# Created on March 25, 2013
#


# source files.
BASE_SRC = cel_timer.c cel_domain_3d.c cel_vector.c  cel_sort.c cel_sparse_matrix_coo.c cel_sparse_matrix_csr.c cel_discretisation_poisson_3d.c cel_discretisation.c 
BASE_SRC += cel_sparce_matrix_coo_interface.c
#output library 
BASE_OUT = $(LIB_DIR)/libcelbase.a

#directories

export BASE_SRC_DIR = $(PROJECT_DIR)/src/base
export BASE_INC_DIR = $(ROOT_INC_DIR)/include/base
export BASE_OBJ_DIR = $(ROOT_OBJ_DIR)/include/base


# include directories
C_INCLUDES_TMP = $(BASE_INC_DIR)
export C_INCLUDES+=$(foreach d, $(C_INCLUDES_TMP), -I$d)

BASE_SOURCES=$(foreach file,$(BASE_SRC),$(BASE_SRC_DIR)/$(file))
BASE_OBJECTS=$(foreach file,$(BASE_SRC),$(BASE_OBJ_DIR)/$(basename $(file)).o)



$(BASE_OBJ_DIR)/%.o:	$(BASE_SRC_DIR)/%.c $(BASE_OBJ_DIR)
	$(CC) $(C_INCLUDES) $(CFLAGS) -c $< -o $@


task_list= $(BASE_OBJECTS) 
ifeq ($(COMPILER),CRAY)
task_list+= cray_timer 
ADD_TIMER_LIB += $(LIB_DIR)/cel_timer_cray.o
endif
task_list+=$(LIB_DIR)




base: $(BASE_OUT)

$(BASE_OUT): $(task_list) 
	ar rcs $(BASE_OUT) $(BASE_OBJECTS) $(ADD_TIMER_LIB)


cray_timer:  $(LIB_DIR)
	gcc -O2 -c  $(BASE_SRC_DIR)/cel_timer_cray.c -o cel_timer_cray.o
	mv cel_timer_cray.o $(LIB_DIR)/

$(BASE_OBJ_DIR) :
	mkdir $(BASE_OBJ_DIR) -p

clean_base:
	rm -rf $(BASE_OBJECTS) $(BASE_OUT)

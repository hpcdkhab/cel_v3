#
# File:   Makefile
# Project CRESTA (see details on https://www.cresta-project.eu)
# Send you email about the bugs to
# The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
# Dmitry Khabi,  (email: khabi@hlrs.de)
# Created on March 25, 2013


# source files
#MPI_OMP_TEST_SRC = cel_cg.f90 
MPI_OMP_TEST_SRC = cel_elm.f90 
#TEST_SRC = cel_mv.f90
#executable file
MPI_OMP_OUT = cel_cg
#OUT = cel_mv
MPI_OMP_LIBS = 
#linked intern libraries
MPI_OMP_LIBS = $(LIB_DIR)/libcelfelmsolver.a
MPI_OMP_LIBS += $(LIB_DIR)/libcelcgalg.a
MPI_OMP_LIBS += $(LIB_DIR)/libcelfop.a
MPI_OMP_LIBS += $(LIB_DIR)/libcelfprec.a
MPI_OMP_LIBS += $(LIB_DIR)/libcelfspmtdistr.a
MPI_OMP_LIBS += $(LIB_DIR)/libcelfspmtconverter.a
MPI_OMP_LIBS += $(LIB_DIR)/libcelfspmt.a
MPI_OMP_LIBS += $(LIB_DIR)/libcelfcomm.a
MPI_OMP_LIBS += $(LIB_DIR)/libcelfbasedistr.a
MPI_OMP_LIBS += $(LIB_DIR)/libcelfbase.a
MPI_OMP_LIBS += $(LIB_DIR)/libcelbase.a

#ext. libraries (cresta wp4 task 1) for the import (IO) of the matrix in FPS format
#MPI_OMP_LIBS += $(TASK4_LIB_DIR)/libwp4ld.a
#MPI_OMP_LIBS += $(TASK4_LIB_DIR)/libwp4read.a
#MPI_OMP_LIBS += $(TASK4_LIB_DIR)/libwp4pa.a
#MPI_OMP_LIBS += $(TASK4_LIB_DIR)/libwp4io.a
#ext. module
MPI_OMP_LIBS += -lm -lrt
ifeq ($(CRAY),YES)
#MPI_OMP_LIBS +=-ldmapp
endif
ifeq ($(PERF_TOOL),CRAYPAT)
#MPI_OMP_LIBS +=-lcraypat
endif
MPI_OMP_LIBS += $(MKL_LIB)

#local directories
export MPI_OMP_TEST_SRC_DIR = ./src
export MPI_OMP_TEST_MOD_DIR = ./mod
export MPI_OMP_TEST_OBJ_DIR = ./obj

.PHONY: directories

# include directories
MPI_OMP_F_INCLUDES_TMP = $(MPI_OMP_TEST_MOD_DIR)  $(FPREC_MOD_DIR)  $(FSPMTDISTR_MOD_DIR) $(FBASEDISTR_MOD_DIR) $(FSPMT_MOD_DIR) $(FBASE_MOD_DIR) $(FCOMM_MOD_DIR) $(FSPMTCONVERTER_MOD_DIR) $(FSPMTDISTR_MOD_DIR) $(FMKL_INC_DIR) $(TASK4_MOD_DIR) 
MPI_OMP_F_INCLUDES_TMP += $(FOMP_MOD_DIR) $(FOP_MOD_DIR) $(FCGALG_MOD_DIR) $(FELMSOLVER_MOD_DIR)
export MPI_OMP_F_INCLUDES+=$(foreach d, $(MPI_OMP_F_INCLUDES_TMP), -I$d)
export MPI_OMP_F_INCLUDES+=$(F_INCLUDES)



MPI_OMP_TEST_SOURCES=$(foreach file,$(MPI_OMP_TEST_SRC),$(MPI_OMP_TEST_SRC_DIR)/$(file))
MPI_OMP_TEST_OBJECTS=$(foreach file,$(MPI_OMP_TEST_SRC),$(MPI_OMP_TEST_OBJ_DIR)/$(basename $(file)).o)


all: clean directories $(MPI_OMP_OUT) 

cel_cg: $(MPI_OMP_OUT)

$(MPI_OMP_TEST_OBJ_DIR)/%.o: $(MPI_OMP_TEST_SRC_DIR)/%.f90  
	$(MPI_F90) $(MPI_OMP_F_INCLUDES) $(FFLAGS) $(F90_MODULE_DIR)  $(MPI_OMP_TEST_MOD_DIR) -c $< -o $@

$(MPI_OMP_OUT): directories  $(MPI_OMP_TEST_OBJECTS) 
	$(MPI_F90) $(FLFLAGS) -o $(MPI_OMP_OUT) $(MPI_OMP_TEST_OBJECTS)  $(MPI_OMP_LIBS)

directories: $(MPI_OMP_TEST_OBJ_DIR)

$(MPI_OMP_TEST_OBJ_DIR):
	mkdir $(MPI_OMP_TEST_OBJ_DIR) $(MPI_OMP_TEST_MOD_DIR) -p

clean: clean_cel_cg

clean_cel_cg:
	rm -rf $(MPI_OMP_TEST_OBJ_DIR) $(MPI_OMP_OUT)



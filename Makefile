#########################
# MAKEFILE FOR arm_test #
#########################

INTELROOT=/opt/intel/composerxe
MKLROOT=$(INTELROOT)/mkl
ARMAROOT=$(HOME)/libs/armadillo
ARPACKROOT=$(HOME)/libs/arpack-ng-3.1.3
CPP=$(INTELROOT)/bin/icpc
CPUTYPE=native
#CPUTYPE=corei7-avx
#CPUTYPE=core2

#FF=ifort                                                                                                                                                    
#FFLAGS = -O3 -m64 -openmp
CPPFLAGS = -m64 -std=c++11 -Wall -DINTEL -DSYMS -wd488 -wd381
#CPPFLAGS = -m64 -openmp -Wall -DINTEL -march=core2
CFLAGS_DEBUG = $(CPPFLAGS) -g -DVERBOSE
CFLAGS_RELEASE = $(CPPFLAGS) -O3 -DNDEBUG -openmp #-align -march=$(CPUTYPE)

INCLUDE = -I$(ARMAROOT)/include -I$(MKLROOT)/include #-I$(HOME)/HDF5/hdf5/include -I$(HOME)/HDF5/SZIP/include #-I$(INTELDIR)/include
#LINCLUDE = -L $(MKLROOT)/lib/intel64
LINCLUDE = -L$(MKLROOT)/lib/intel64 #-L$(HOME)/lib

# static linking
#LDFLAGS_DEBUG = -wd10237  -static-intel -Wl,--start-group -larpackng_313_debug -lifcore -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread -lm 
LDFLAGS_DEBUG = -wd10237  -static-intel  -Wl,--start-group ${ARPACKROOT}/libarpack_debug.a ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -lgfortran
#LDFLAGS_RELEASE = -wd10237 -static-intel -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread -lm 
LDFLAGS_RELEASE = -wd10237 -static-intel  -Wl,--start-group ${ARPACKROOT}/libarpack_release.a ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -lpthread -lm -liomp5 -lgfortran
# dynamic linking
#LDFLAGS = 

OBJECTS = ITEBD_main.o
EXENAME = ITEBD_mkl

OBJDIR_DEBUG = obj/DebugMKL
OUT_DEBUG = bin/Debug/$(EXENAME)
OBJ_DEBUG = $(OBJECTS:%.o=$(OBJDIR_DEBUG)/%.o)

OBJDIR_RELEASE = obj/ReleaseMKL
OUT_RELEASE = bin/Release/$(EXENAME)
OBJ_RELEASE = $(OBJECTS:%.o=$(OBJDIR_RELEASE)/%.o)

HELP_DIR = ../../helpers
HELP_OBJ_DIR_DEBUG = $(OBJDIR_DEBUG)/helpers
HELP_OBJ_DIR_RELEASE = $(OBJDIR_RELEASE)/helpers

HELP_HEAD = $(HELP_DIR)/helpers.hpp $(HELP_DIR)/arma_typedefs.h $(HELP_DIR)/tictoc.hpp $(HELP_DIR)/BlockObj.hpp $(HELP_DIR)/IKey.hpp $(HELP_DIR)/MPSBlockMat.hpp $(HELP_DIR)/MPSBlockUtilities.hpp $(HELP_DIR)/EigsBlockFunctions.hpp
HELPERS = parser.o eigs.o
HELPERS_SRC = $(HELPERS:%.o=%.h)
HELP_OBJ_DEBUG = $(HELPERS:%.o=$(HELP_OBJ_DIR_DEBUG)/%.o)
HELP_OBJ_RELEASE = $(HELPERS:%.o=$(HELP_OBJ_DIR_RELEASE)/%.o)
#-----------------------------------------------------------------------------------
all: debug release
#-----------------------------------------------------------------------------------
test:
	@echo '$(HELP_OBJ_DEBUG)'
	@echo '$(HELP_OBJ_RELEASE)'
	@echo '$(LINCLUDE)'
#-----------------------------------------------------------------------------------
debug: before_debug out_debug

before_debug:
	@echo 'Checking DIRS for debug:'
	test -d bin/Debug || mkdir -p bin/Debug
	test -d $(OBJDIR_DEBUG) || mkdir -p $(OBJDIR_DEBUG)
	test -d $(HELP_OBJ_DIR_DEBUG) || mkdir -p $(HELP_OBJ_DIR_DEBUG)
	@echo 'done'
	@echo ' '

$(OBJDIR_DEBUG)/ITEBD_main.o: ITEBD_main.cpp $(HELP_HEAD)
	@echo 'Building file: $@'
	$(CPP) $(CFLAGS_DEBUG) $(INCLUDE) -c $< -o $@
	@echo 'Finished building: $@'
	@echo ' '

$(HELP_OBJ_DIR_DEBUG)/%.o: $(HELP_DIR)/%.cpp $(HELP_DIR)/%.h $(HELP_HEAD)
	@echo 'Building file: $@'
	$(CPP) $(CFLAGS_DEBUG) $(INCLUDE) -c $< -o $@
	@echo 'Finished building: $@'
	@echo ' '

out_debug: $(OBJ_DEBUG) $(HELP_OBJ_DEBUG)
	@echo 'Linking: $(OUT_DEBUG)'
	$(CPP) -o $(OUT_DEBUG) $^ $(LINCLUDE) $(LDFLAGS_DEBUG)
	@echo 'Finished linking $(OUT_DEBUG)'
	@echo ' '
#-----------------------------------------------------------------------------------
release: before_release out_release

before_release:
	@echo 'Checking DIRS for release:'
	test -d bin/Release || mkdir -p bin/Release
	test -d $(OBJDIR_RELEASE) || mkdir -p $(OBJDIR_RELEASE)
	test -d $(HELP_OBJ_DIR_RELEASE) || mkdir -p $(HELP_OBJ_DIR_RELEASE)
	@echo 'done'
	@echo ' '

$(OBJDIR_RELEASE)/ITEBD_main.o: ITEBD_main.cpp $(HELP_HEAD)
	@echo 'Building file: $@'
	$(CPP) $(CFLAGS_RELEASE) $(INCLUDE) -c $< -o $@
	@echo 'Finished building: $@'
	@echo ' '

$(HELP_OBJ_DIR_RELEASE)/%.o: $(HELP_DIR)/%.cpp $(HELP_DIR)/%.h $(HELP_HEAD)
	@echo 'Building file: $@'
	$(CPP) $(CFLAGS_RELEASE) $(INCLUDE) -c $< -o $@
	@echo 'Finished building: $@'
	@echo ' '

out_release: $(OBJ_RELEASE) $(HELP_OBJ_RELEASE)
	@echo 'Linking: $(OUT_RELEASE)'
	$(CPP) -o $(OUT_RELEASE) $^ $(LINCLUDE) $(LDFLAGS_RELEASE)
	@echo 'Finished linking $(OUT_RELEASE)'
	@echo ' '

#-----------------------------------------------------------------------------------
clean: cleanDebug cleanRelease

cleanDebug: 
	rm -f $(HELP_OBJ_DEBUG) $(OBJ_DEBUG) $(OUT_DEBUG)
	rm -rf bin/Debug
	rm -rf $(OBJDIR_DEBUG)

cleanRelease: 
	rm -f $(HELP_OBJ_RELEASE) $(OBJ_RELEASE) $(OUT_RELEASE)
	rm -rf bin/Release
	rm -rf $(OBJDIR_RELEASE)

.PHONY: before_debug cleanDebug before_release cleanRelease test

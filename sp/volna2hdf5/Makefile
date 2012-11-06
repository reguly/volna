#
# The following environment variables should be predefined:
#
# OP2_INSTALL_PATH
# HDF5_INSTALL_PATH
# VOLNA_INSTALL_PATH
# MPI_INSTALL_PATH
#

#
# set paths for header files and libraries
#

VOLNA_INSTALL_PATH ?= ../../volna_init

OP2_INC = -I$(OP2_INSTALL_PATH)/c/include
OP2_LIB = -L$(OP2_INSTALL_PATH)/c/lib

VOLNA_INC = -I$(VOLNA_INSTALL_PATH)
VOLNA_SRC = $(VOLNA_INSTALL_PATH)

HDF5_INC = -I$(HDF5_INSTALL_PATH)/include
HDF5_LIB = -L$(HDF5_INSTALL_PATH)/lib -lhdf5_hl -lhdf5 -lz

MPI_INC = -I$(MPI_INSTALL_PATH)/include

ifeq ($(OP2_COMPILER),gnu)
	CPP = g++
	CPPFLAGS = -O3 -msse4.2 -fPIC -DUNIX -Wall -DOP2 -arch x86_64 -fopenmp
	#CPPFLAGS = -g -O0 -msse4.2 -fPIC -DUNIX -Wall -DOP2 -arch x86_64 -fopenmp
	MPICPP = mpic++
	MPIFLAGS = $(CCFLAGS)
else
ifeq ($(OP2_COMPILER),intel)
	CPP = icpc
	#CPPFLAGS = -g -O0 -vec-report -xSSE4.2 -parallel -openmp -openmp-report2
	CPPFLAGS = -O3 -vec-report -xSSE4.2 -parallel -openmp -openmp-report2
	MPICPP = mpiCC
	MPIFLAGS = $(CPPFLAGS)
else
print:
		@echo "unrecognized value for OP2_COMPILER"
endif
endif


#
# master to make all versions
#


VOLNA_SOURCES = \
	$(VOLNA_SRC)/values.cpp\
	$(VOLNA_SRC)/event.cpp
#  $(VOLNA_SRC)/geom.cpp\

VOLNA_OBJECTS = $(VOLNA_SOURCES:.cpp=.o)

all: clean volna2hdf5

volna2hdf5: volna2hdf5.cpp $(VOLNA_OBJECTS) Makefile
	$(MPICPP) $(CPPFLAGS) volna2hdf5.cpp $(VOLNA_OBJECTS) -I$(VOLNA_INSTALL_PATH) \
            -I$(VOLNA_INSTALL_PATH)/external $(VOLNA_INC) $(HDF5_INC) $(OP2_INC) \
            $(OP2_LIB) $(HDF5_LIB) -lop2_seq -lop2_hdf5 -o volna2hdf5

#
# cleanup
#
clean:
		rm -f volna2hdf5 *.o $(VOLNA_OBJECTS)

%.o: %.cpp %.hpp
	$(CPP) $(CPPFLAGS) -I$(VOLNA_INSTALL_PATH) -I$(VOLNA_INSTALL_PATH)/external -c -o $@ $<

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

VOLNA_INSTALL_PATH = .

LIBMESH_INC = -I$(LIBMESH_INSTALL_PATH)/include
LIBMESH_LIB = -L$(LIBMESH_INSTALL_PATH)/lib

OP2_INC = -I$(OP2_INSTALL_PATH)/c/include
OP2_LIB = -L$(OP2_INSTALL_PATH)/c/lib

VOLNA_INC = -I$(VOLNA_INSTALL_PATH)
VOLNA_SRC = $(VOLNA_INSTALL_PATH)

HDF5_INC = -I$(HDF5_INSTALL_PATH)/include
HDF5_LIB = -L$(HDF5_INSTALL_PATH)/lib -lhdf5_hl -lhdf5 -lz

MPI_INC = -I$(MPI_INSTALL_PATH)/include

ifeq ($(OP2_COMPILER),gnu)
	CPP = g++
#	CPPFLAGS = -g0 -msse3 -fPIC -DUNIX -Wall
	CPPFLAGS = -g -O0 -msse4.2 -fPIC -DUNIX -Wall -DOP2
	#CPPFLAGS = -g -fPIC -DUNIX -Wall
	OMPFLAGS = -fopenmp
	MPICPP = mpic++
	MPIFLAGS = $(CCFLAGS)
else
ifeq ($(OP2_COMPILER),intel)
	CPP = icpc
	CCFLAGS = -g -O0 -vec-report -xSSE4.2 -parallel
	#CCFLAGS = -g -vec-report -xSSE4.2 -parallel
	CPPFLAGS = $(CCFLAGS)
	OMPFLAGS = -openmp -openmp-report2
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

VOLNA_OBJECTS = $(VOLNA_SOURCES:.cpp=.o)

all: clean $(VOLNA_OBJECTS)

#
# cleanup
#

clean:
		rm -f $(VOLNA_OBJECTS)

%.o: %.cpp %.hpp
	$(CPP) $(CPPFLAGS) -I$(VOLNA_INSTALL_PATH) -c -o $@ $<

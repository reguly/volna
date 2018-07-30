#
# The following environment variables should be predefined:
#
# OP2_INSTALL_PATH =
# HDF5_INSTALL_PATH =
# MPI_INSTALL_PATH =
# PARMETIS_INSTALL_PATH =
# PTSCOTCH_INSTALL_PATH =
# CUDA_INSTALL_PATH =
# VOLNA_INSTALL_PATH =
# OP2_COMPILER =

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

ifdef PTSCOTCH_INSTALL_PATH
	PTSCOTCH_INC 	= -I$(PTSCOTCH_INSTALL_PATH)/include -DHAVE_PTSCOTCH
	PTSCOTCH_LIB 	= -L$(PTSCOTCH_INSTALL_PATH)/lib/ -lptscotch \
	-L$(PTSCOTCH_INSTALL_PATH)/lib/ -lptscotcherr
endif

ifeq ($(OP2_COMPILER),gnu)
	CPP = g++
	#CPPFLAGS = -O2 -msse4.2 -fPIC -DUNIX -Wall -DOP2 -fopenmp -fpermissive
	#CPPFLAGS = -g -O0 -msse4.2 -fPIC -DUNIX -Wall -DOP2 -arch x86_64 -fopenmp
	CPPFLAGS = -O2  -fPIC -DUNIX -Wall -DOP2 -fopenmp -DEIGEN_DONT_VECTORIZE -DEIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT
	MPICPP = $(MPI_INSTALL_PATH)/bin/mpic++
	MPIFLAGS = $(CCFLAGS)
else
ifeq ($(OP2_COMPILER),intel)
	CPP = icpc
	CPPFLAGS = -O2 -g -parallel -qopenmp  -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX
	#CPPFLAGS = -g -O0 -vec-report -xSSE4.2 -parallel -openmp -openmp-report2
	MPICPP = $(MPI_INSTALL_PATH)/bin/mpiCC
	MPIFLAGS = $(CPPFLAGS)
else
print:
		@echo "unrecognized value for OP2_COMPILER"
endif
endif


VOLNA_SOURCES = \
	$(VOLNA_SRC)/values.cpp\
	$(VOLNA_SRC)/event.cpp

VOLNA_OBJECTS = $(VOLNA_SOURCES:.cpp=.o)

all: clean volna2hdf5

volna2hdf5: volna2hdf5.cpp ../volna_util.cpp $(VOLNA_OBJECTS) Makefile
ifdef MPI_INSTALL_PATH
	$(MPICPP) $(CPPFLAGS) volna2hdf5.cpp ../volna_util.cpp $(VOLNA_OBJECTS) $(PTSCOTCH_INC) $(PTSCOTCH_LIB) -lscotch -lscotcherr -I$(VOLNA_INSTALL_PATH) \
            -I$(VOLNA_INSTALL_PATH)/external $(VOLNA_INC) $(HDF5_INC) $(OP2_INC) \
            $(OP2_LIB) $(HDF5_LIB) -lop2_seq -lop2_hdf5 -o volna2hdf5
else
	$(CPP) $(CPPFLAGS) volna2hdf5.cpp ../volna_util.cpp $(VOLNA_OBJECTS) $(PTSCOTCH_INC) $(PTSCOTCH_LIB) -lscotch -I$(VOLNA_INSTALL_PATH) \
            -I$(VOLNA_INSTALL_PATH)/external $(VOLNA_INC) $(HDF5_INC) $(OP2_INC) \
            $(OP2_LIB) $(HDF5_LIB) -lop2_seq -lop2_hdf5 -o volna2hdf5
endif

#
# cleanup
#
clean:
		rm -f volna2hdf5 *.o $(VOLNA_OBJECTS)

%.o: %.cpp %.hpp
	$(CPP) $(CPPFLAGS) -I$(VOLNA_INSTALL_PATH) -I$(VOLNA_INSTALL_PATH)/external -c -o $@ $<

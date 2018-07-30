#
# The following environment variables should be predefined:
OP2_INSTALL_PATH= /usr/local/OP2-Common/op2
HDF5_INSTALL_PATH= /usr/local/hdf5-1.8.15/hdf5
MPI_INSTALL_PATH= /usr/lib64/openmpi
#

#
# set paths for header files and libraries
#
OP2_COMPILER= gnu

OP2_INC = -I$(OP2_INSTALL_PATH)/c/include
OP2_LIB = -L$(OP2_INSTALL_PATH)/c/lib

HDF5_INC = -I$(HDF5_INSTALL_PATH)/include
HDF5_LIB = -L$(HDF5_INSTALL_PATH)/lib 

ifeq ($(OP2_COMPILER),gnu)
	CPP = g++
	CPPFLAGS = -O2 -msse4.2 -fPIC -DUNIX -Wall -DOP2 -fopenmp
	#CPPFLAGS = -g -O0 -msse4.2 -fPIC -DUNIX -Wall -DOP2 -arch x86_64 -fopenmp
	MPICPP = $(MPI_INSTALL_PATH)/bin/mpic++
	MPIFLAGS = $(CCFLAGS)
else
ifeq ($(OP2_COMPILER),intel)
	CPP = icpc
	CPPFLAGS = -O2 -vec-report -xSSE4.2 -parallel -openmp -openmp-report2 -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX
	#CPPFLAGS = -O2 -vec-report -xSSE4.2 -parallel -openmp -openmp-report2
	MPICPP = $(MPI_INSTALL_PATH)/bin/mpiCC
	MPIFLAGS = $(CPPFLAGS)
else
print:
		@echo "unrecognized value for OP2_COMPILER"
endif
endif


#
# master to make all versions
#
all: clean hdf52vtk

hdf52vtk: hdf52vtk.cpp ../volna_writeVTK.cpp Makefile
ifdef MPI_INSTALL_PATH
	$(MPICPP) $(CPPFLAGS) $(HDF5_LIB) $(OP2_LIB) ../volna_writeVTK.cpp hdf52vtk.cpp -lop2_hdf5 -lhdf5 -lz -lop2_seq -lhdf5_hl $(HDF5_INC) $(OP2_INC) \
            -o hdf52vtk
else
	$(CPP) $(CPPFLAGS) $(HDF5_LIB) $(OP2_LIB) ../volna_writeVTK.cpp hdf52vtk.cpp -lop2_hdf5 -lop2_seq $(HDF5_INC) $(OP2_INC) \
            -o hdf52vtk 
endif

#
# cleanup
#
clean:
		rm -f hdf52vtk *.o 


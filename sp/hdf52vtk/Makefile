#
# The following environment variables should be predefined:
#
# OP2_INSTALL_PATH
# HDF5_INSTALL_PATH
# MPI_INSTALL_PATH - (optional) if use of MPI is considered 
#

#
# set paths for header files and libraries
#

OP2_INC = -I$(OP2_INSTALL_PATH)/c/include
OP2_LIB = -L$(OP2_INSTALL_PATH)/c/lib


HDF5_INC = -I$(HDF5_INSTALL_PATH)/include
HDF5_LIB = -L$(HDF5_INSTALL_PATH)/lib 

ifeq ($(OP2_COMPILER),gnu)
	CPP = g++
	CPPFLAGS = -O2 -msse4.2 -fPIC -DUNIX -Wall -DOP2 -arch x86_64 -fopenmp
	#CPPFLAGS = -g -O0 -msse4.2 -fPIC -DUNIX -Wall -DOP2 -arch x86_64 -fopenmp
	MPICPP = mpic++
	MPIFLAGS = $(CCFLAGS)
else
ifeq ($(OP2_COMPILER),intel)
	CPP = icpc
	CPPFLAGS = -g -O0 -vec-report -xSSE4.2 -parallel -openmp -openmp-report2
	#CPPFLAGS = -O2 -vec-report -xSSE4.2 -parallel -openmp -openmp-report2
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

all: clean hdf52vtk

hdf52vtk: hdf52vtk.cpp ../volna_writeVTK.cpp Makefile
ifdef MPI_INSTALL_PATH
	$(MPICPP) $(CPPFLAGS) $(HDF5_LIB) $(OP2_LIB) ../volna_writeVTK.cpp hdf52vtk.cpp -lop2_hdf5 -lhdf5 -lz -lop2_seq -lhdf5_hl $(HDF5_INC) $(OP2_INC) \
            -o hdf52vtk
else
	$(CPP) $(CPPFLAGS) $(HDF5_LIB) $(OP2_LIB) ../volna_writeVTK.cpp hdf52vtk.cpp -lop2_hdf5 -lop2_seq $(HDF5_INC) $(OP2_INC) \
            -o hdf52vtk 
#	$(CPPFLAGS) hdf52vtk.cpp ../volna_output.cpp $(HDF5_INC) $(OP2_INC) \
#            $(OP2_LIB) $(HDF5_LIB) -lop2_seq -lop2_hdf5 -o hdf52vtk
endif

#
# cleanup
#
clean:
		rm -f hdf52vtk *.o 


Volna- OP2
=====

This is an OP2 port of the original Volna code.

## Installation
You need to install OP2 (source in [this](https://github.com/OP2/OP2-Common) GitHub repository)
 * You will need HDF5 (preferably compiled with MPI support)
 * You may not need the partitioners PT-Scotch/ParMetis if you do not plan to run in a distributed environment

Check out the Volna-OP2 code from [this](https://github.com/reguly/Volna) GitHub repository. The code comprises of two parts, volna2hdf5 which takes configuration files from the original Volna code and dumps all the necessary information to a h5 file which will be used by the second application volna-OP2.
 * Type 'make' in sp/volna2hdf5 - note that some warnings will show because of dependencies, you can safely ignore these
 * Type make in sp to build all versions of the volna-op2 application
  * make volna_seq builds the sequential (single-threaded) version
  * make volna_openmp builds the single node multi-threaded version
  * make volna_cuda builds the single node GPU version
  * make volna_mpi builds the sequential MPI version with a single thread executing on each MPI process
  * make volna_mpi_openmp builds the MPI+OpenMP version
  * make volna_mpi_cuda builds the MPI+CUDA version
		
## Use
To use volna-OP2 with the *.vln configuration files, first you have to use volna2hdf5, e.g.
 * ./volna2hdf5 gaussian_landslide.vln which will output a gaussian_landslide.h5 file
Afterwards, call volna-op2 with the above input file, e.g.:
 * ./volna_openmp gaussian_landslide.h5
 * when using the CUDA version we suggest adding "OP_PART_SIZE=128 OP_BLOCK_SIZE=128" to the execution line

## Output
The files written by volna-OP2 may be different from the original VOLNA code, due to performance optimisations.
 * OutputLocation events are bundled together, if they all use the same timings, into a gauges.h5 file, which is compressed to save space. To read the file, you can use any hdf5 tool, there are two datasets written: /dims which has 2 integer fields, the first indicating the size in the contiguous direction (number of OutputLocation events + 1 for timestaps) and the second in the strided dimension (number of timestamps). Data is stored under /gauges, a 1D float array of size dims[0]*dims[1]. To get the timestamp at iteration T and the value for gauge N (indexed 1...dims[0]-1), access gauges[T*dims[0]] for the timestamp and gauges[T*dims[0]+N] - T has to be less than dims[1]. A python script is attached (read_gauges.py) that shows an example of this.

## Recommendations, restrictions
Some restriction, constantly updated as they are fixed:
 * Currently Volna does not work well with MPI, and there is no support for distributed file output (e.g. OutputSimulation or OutputLocation into files).
 * Refrain from using OutputSimulation too often because (compared to the actual simulation) it may take a lot of time.
 * When using large h5 files, try compressing them, e.g. h5repack -i gaussian.h5 -o gaussian_compressed.h5 -f GZIP=9

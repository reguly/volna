HDF52VTK - tool to convert results to VTK

Syntax: hdf52vtk geometry_filename.h5 sim_filename.h5 [0,1] [output.vtk]

The outputs of Volna OP2 are HDF5 files with *.h5 extension. These files contain the 
snapshots of simulation variables. Using this tool the raw OP2 data can be converted
to VTK files. One has to specify at least to files: one containing the geometry data
and one containing the simulation data. For consecutive simulation outputs (*.h5) 
the geometry data is the same.

geometry_filename.h5 - Contains the geometric data that is used by the simulation.
                       Created by: ./volna2hdf5 script_filename.vln
sim_result.h5        - Contains the simulation results.
                       Created by: ./volna geometry_filename.h5
0,1                  - [optional] 0(default) - write data to ASCII VTK; 1 - write data to Binary VTK.
output.vtk           - [optional] The name of the output VTK file. By default output.vtk

The data of the two files (geometry_filename and sim_result) are combined and are put into *.vtk files\

Examples:
  $ hdf52vtk cata2.h5 sim_result.h5
  $ hdf52vtk cata2.h5 sim_result.h5 0 output.vtk

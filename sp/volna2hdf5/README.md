volna2hdf5 - save Volna data to HDF5 file
==========

    Synthax: volna2hdf5 <filename.vln> [no-reorder]

Transfers Volna specific data to OP2 HDF5 file. The Volna config file with *.vln extension has to be specified. The tool uses the given config file to produce an HDF5 file that contains all the data necessary to run the OP2 port of Volna. The produced HDF5 file has the same file name with *.h5 extension. During the HDF5 file creation maps and data elements (cells, nodes, edges and associated data) are reordered to make the run-time more efficient. If no reordering is needed (e.g. in case of comparison with original Volna results) the no-reorder option can be used. 

Examples:
    $ ./volna2hdf5 cata2.vln 
    $ ./volna2hdf5 cata2.vln no-reorder

If necessary, the HDF5 file can be viewed using h5dump:  

    $ h5dump file.h5 > log && vim log


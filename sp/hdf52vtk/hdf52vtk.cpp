#include"../volna_common.h" // macro definitions
#include<math.h>
#include "../volna_writeVTK.h"

//
// Checks if error occured during hdf5 process and prints error message
//
void __check_hdf5_error(herr_t err, const char *file, const int line){
  if (err < 0) {
    printf("%s(%i) : OP2_HDF5_error() Runtime API error %d.\n", file,
        line, (int) err);
    exit(-1);
  }
}

void print_info() {
  printf("Wrong parameters! Please specify the OP2 HDF5 data filenames. \n"
    "Syntax: hdf52vtk geometry_filename.h5 sim_result.h5 [0,1] output.vtk\n\n"
    "  geometry_filename.h5 - Contains the geometric data that is used by the simulation.\n"
    "                         Created by: ./volna2hdf5 script_filename.vln\n"
    "  sim_result.h5        - Contains the simulation results.\n"
    "                         Created by: ./volna geometry_filename.h5\n\n"
    "  0,1                  - 0(default) - write data to ASCII VTK; 1 - write data to Binary VTK.\n"
    "  output.vtk           - The name of the output VTK file. By default output.vtk\n"
    "The data of the two files (geometry_filename and sim_result) are combined and are put into *.vtk files\n");
//    printf("Wrong parameters! Please specify the OP2 HDF5 data filenames. \n"
//      "Syntax: hdf52vtk geometry_filename.h5 sim_result.h5 [sim_diff.h5] \n\n"
//      "  geometry_filename.h5 - contains the geometric data that is used by the simulation.\n"
//      "                         Created by: ./volna2hdf5 script_filename.vln\n"
//      "  sim_result.h5        - contains the simulation results.\n"
//      "                         Created by: ./volna geometry_filename.h5\n\n"
//      "  sim_diff.h5          - [Optional] Used to calculate differences along the 4 dimesions."
//      "                         Contains simulation results.\n"
//      "                         Created by: ./volna geometry_filename.h5\n\n"
//      "The data of the two files (geometry_filename and sim_result) are combined and are put into *.vtk files\n");
}

int main(int argc, char **argv) {
  if (argc < 3) {
    print_info();
    exit(-1);
  }
  op_init(argc, argv, 2);
  hid_t file_geo;
  hid_t file_sim;
  const char *filename_geo_h5 = argv[1];
  const char *filename_sim_h5 = argv[2];
  file_geo = H5Fopen(filename_geo_h5, H5F_ACC_RDONLY, H5P_DEFAULT);
  file_sim = H5Fopen(filename_sim_h5, H5F_ACC_RDONLY, H5P_DEFAULT);

  // Define OP2 sets - Read mesh and geometry data from HDF5
  op_set nodes = op_decl_set_hdf5(filename_geo_h5, "nodes");
  op_set cells = op_decl_set_hdf5(filename_geo_h5, "cells");

  // Define OP2 set maps
  op_map cellsToNodes = op_decl_map_hdf5(cells, nodes, N_NODESPERCELL,
  		filename_geo_h5,
  		"cellsToNodes");

  // Define OP2 set dats
  op_dat nodeCoords = op_decl_dat_hdf5(nodes, MESH_DIM, "float",
  		filename_geo_h5,
  		"nodeCoords");
  op_dat physical_vars = op_decl_dat_hdf5(cells, N_STATEVAR+1, "float",
  		filename_sim_h5,
  		"physical_vars");

  int nnode = nodeCoords->set->size;
  int ncell = cellsToNodes->from->size;

	// Set output filename
	char filename[255];
	if(argc == 5) {
	  sprintf(filename, "%s", argv[4]);
	} else {
	  sprintf(filename, "%s", filename_sim_h5);
	  const char* substituteIndexPattern = ".h5";
	  char* pos;
	  pos = strstr(filename, substituteIndexPattern);
	  char substituteIndex[255];
	  sprintf(substituteIndex, ".vtk");
	  strcpy(pos, substituteIndex);
	}
	// Choose file type
	if(argc == 4) {
	  switch( atoi(argv[3]) ) {
	  case 0:
	    WriteMeshToVTKAscii(filename, nodeCoords, nnode, cellsToNodes, ncell, physical_vars);
	    break;
	  case 1:
	    WriteMeshToVTKBinary(filename, nodeCoords, nnode, cellsToNodes, ncell, physical_vars);
	    break;
	  default:
	    print_info();
	    exit(-1);
	    break;
	  }
	}
	else {
	  WriteMeshToVTKAscii(filename, nodeCoords, nnode, cellsToNodes, ncell, physical_vars);
	}

  check_hdf5_error( H5Fclose(file_geo) );
  check_hdf5_error( H5Fclose(file_sim) );
  op_exit();
  return 0;
}

/*Copyright 2018, Frederic Dias, Serge Guillas, Istvan Reguly

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "volna_common.h"
#include "volna_util.h"
#include "EvolveValuesRK2_1.h"
#include "EvolveValuesRK2_2.h"
#include "simulation_1.h"
#include "limits.h"
#include "Friction_manning.h"
#include "zero_bathy.h"
//
// Sequential OP2 function declarations
//
#include "op_seq.h"

//these are not const, we just don't want to pass them around
LocationData locationData;
double timestamp = 0.0;
int itercount = 0;

// Constants
float CFL, g, EPS;
bool new_format = true;

// Store maximum elevation and speed in global variable, for the sake of max search
op_dat currentMaxElevation = NULL;
op_dat currentMaxSpeed = NULL;
op_dat physical_vars = NULL;
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
  op_printf("\nPlease specify the OP2 HDF5 data file "
        "name, which was created with volna2hdf5 tool, e.g. volna2hdf5 script_filename.vln. \n"
        "Use Volna configuration script filename with the *.h5 extension and "
        "specify the output filetype: 0 - HDF5, 1 - VTK ASCII, 2 - VTK Binary \n"
        "e.g. ./volna script_filename.h5 1 \n"
        "or   ./volna script_filename.h5 0 \n"
        "or   mpirun -np 8 ./volna_mpi script_filename.h5 1 \n");
}

int main(int argc, char **argv) {
  op_init(argc, argv, 2);
  if(argc < 3) {
    op_printf("Wrong parameters! \n");
    print_info();
    exit(-1);
  }

  for ( int n = 1; n < argc; n++ )
  {
    if ( strncmp ( argv[n], "old-format", 10 ) == 0 ) {
      new_format = false;
      op_printf("Using old format bathymetry\n");
    }
  }

  const char *filename_h5 = argv[1];
  int writeOption = atoi(argv[2]); // 0 - HDF5, 1 - VTK ASCII, 2 - VTK Binary

  EPS = 1e-6; //machine epsilon, for doubles 1e-11

  hid_t file;
  file = H5Fopen(filename_h5, H5F_ACC_RDONLY, H5P_DEFAULT);
  //Some simulation parameters when using InitGaussianLandslide and InitBore
  GaussianLandslideParams gaussian_landslide_params;
  BoreParams bore_params;
  //Parameters for the rectangualr domain special case
  RectangleDomainParams rect_params;
  //Read the above parameters
  check_hdf5_error( H5LTread_dataset_float(file, "BoreParamsx0", &bore_params.x0) );
  check_hdf5_error( H5LTread_dataset_float(file, "BoreParamsHl", &bore_params.Hl) );
  check_hdf5_error( H5LTread_dataset_float(file, "BoreParamsul", &bore_params.ul) );
  check_hdf5_error( H5LTread_dataset_float(file, "BoreParamsvl", &bore_params.vl) );
  check_hdf5_error( H5LTread_dataset_float(file, "BoreParamsS", &bore_params.S) );
  check_hdf5_error( H5LTread_dataset_float(file, "GaussianLandslideParamsA", &gaussian_landslide_params.A) );
  check_hdf5_error( H5LTread_dataset_float(file, "GaussianLandslideParamsv", &gaussian_landslide_params.v) );
  check_hdf5_error( H5LTread_dataset_float(file, "GaussianLandslideParamslx", &gaussian_landslide_params.lx) );
  check_hdf5_error( H5LTread_dataset_float(file, "GaussianLandslideParamsly", &gaussian_landslide_params.ly) );
  check_hdf5_error( H5LTread_dataset_int(file, "nx", &rect_params.nx) );
  check_hdf5_error( H5LTread_dataset_int(file, "ny", &rect_params.ny) );
  check_hdf5_error( H5LTread_dataset_float(file, "xmin", &rect_params.xmin) );
  check_hdf5_error( H5LTread_dataset_float(file, "xmax", &rect_params.xmax) );
  check_hdf5_error( H5LTread_dataset_float(file, "ymin", &rect_params.ymin) );
  check_hdf5_error( H5LTread_dataset_float(file, "ymax", &rect_params.ymax) );

  int num_events = 0;
  int num_outputLocation = 0;

  check_hdf5_error(H5LTread_dataset_int(file, "numEvents", &num_events));
  std::vector<TimerParams> timers(num_events);
  std::vector<EventParams> events(num_events);

	//Read Event "objects" (Init and Output events) into timers and events
  read_events_hdf5(file, num_events, &timers, &events, &num_outputLocation);

  check_hdf5_error(H5Fclose(file));

  /*
   * Define OP2 sets - Read mesh and geometry data from HDF5
   */
  op_set nodes = op_decl_set_hdf5(filename_h5, "nodes");
  op_set edges = op_decl_set_hdf5(filename_h5, "edges");
  op_set cells = op_decl_set_hdf5(filename_h5, "cells");

  /*
   * Define OP2 set maps
   */
  op_map cellsToCells = op_decl_map_hdf5(cells, cells, N_NODESPERCELL,
      filename_h5,
      "cellsToCells");
  op_map cellsToNodes = op_decl_map_hdf5(cells, nodes, N_NODESPERCELL,
      filename_h5,
      "cellsToNodes");
  op_map edgesToCells = op_decl_map_hdf5(edges, cells, N_CELLSPEREDGE,
      filename_h5,
      "edgesToCells");
  op_map cellsToEdges = op_decl_map_hdf5(cells, edges, N_NODESPERCELL,
      filename_h5,
      "cellsToEdges");

  // When using OutputLocation events we have already computed the cell
  // index of the points so we don't have to locate the cell every time
  op_set outputLocation = NULL;
  op_map outputLocation_map = NULL;
  op_dat outputLocation_dat = NULL;
  if (num_outputLocation) {
    outputLocation = op_decl_set_hdf5(filename_h5, "outputLocation");
    outputLocation_map = op_decl_map_hdf5(outputLocation, cells, 1,
        filename_h5,
        "outputLocation_map");
  }

  /*
   * Define OP2 datasets
   */
  op_dat cellCenters = op_decl_dat_hdf5(cells, MESH_DIM, "float",
                                    filename_h5,
                                    "cellCenters");

  op_dat edgeCenters = op_decl_dat_hdf5(edges, MESH_DIM, "float",
                                    filename_h5,
                                    "edgeCenters");

  op_dat cellVolumes = op_decl_dat_hdf5(cells, 1, "float",
                                    filename_h5,
                                    "cellVolumes");

  op_dat edgeNormals = op_decl_dat_hdf5(edges, MESH_DIM, "float",
                                    filename_h5,
                                    "edgeNormals");

  op_dat edgeLength = op_decl_dat_hdf5(edges, 1, "float",
                                    filename_h5,
                                    "edgeLength");

  op_dat nodeCoords = op_decl_dat_hdf5(nodes, MESH_DIM, "float",
                                      filename_h5,
                                      "nodeCoords");

  op_dat values = op_decl_dat_hdf5(cells, N_STATEVAR, "float",
                                    filename_h5,
                                    "values");
  op_dat isBoundary = op_decl_dat_hdf5(edges, 1, "int",
                                    filename_h5,
                                    "isBoundary");

  //op_dats storing InitBathymetry and InitEta event files
  op_dat temp_initEta         = NULL;
  op_dat temp_initU         = NULL;
  op_dat temp_initV         = NULL;
  op_dat* temp_initBathymetry = NULL;  // Store initBathymtery in an array: there might be more input files for different timesteps
  int n_initBathymetry = 0; // Number of initBathymetry files
  op_set bathy_nodes;
  op_set lifted_cells;
  op_map liftedcellsToBathyNodes;
  op_map liftedcellsToCells;
  op_dat bathy_xy;
  op_dat initial_zb = NULL;


  /*
   * Read constants from HDF5
   */
  op_get_const_hdf5("CFL", 1, "float", (char *) &CFL, filename_h5);

  // Final time: as defined by Volna the end of real-time simulation
  float ftime_tmp, dtmax_tmp;
  op_get_const_hdf5("ftime", 1, "float", (char *) &ftime_tmp, filename_h5);
  op_get_const_hdf5("dtmax", 1, "float", (char *) &dtmax_tmp, filename_h5);
  double ftime = ftime_tmp;
  double dtmax = dtmax_tmp;
  op_get_const_hdf5("g", 1, "float", (char *) &g, filename_h5);

  op_decl_const(1, "float", &CFL);
  op_decl_const(1, "float", &EPS);
  op_decl_const(1, "float", &g);

  //Read InitBathymetry and InitEta event data when they come from files
  for (unsigned int i = 0; i < events.size(); i++) {
      if (!strcmp(events[i].className.c_str(), "InitEta")) {
        if (strcmp(events[i].streamName.c_str(), ""))
          temp_initEta = op_decl_dat_hdf5(cells, 1, "float",
              filename_h5,
              "initEta");
      } else if (!strcmp(events[i].className.c_str(), "InitU")) {
        if (strcmp(events[i].streamName.c_str(), ""))
          temp_initU = op_decl_dat_hdf5(cells, 1, "float",
              filename_h5,
              "initU");
      } else if (!strcmp(events[i].className.c_str(), "InitV")) {
        if (strcmp(events[i].streamName.c_str(), ""))
          temp_initV = op_decl_dat_hdf5(cells, 1, "float",
              filename_h5,
              "initV");
      } else if (!strcmp(events[i].className.c_str(), "InitBathymetry")) {
        if (strcmp(events[i].streamName.c_str(), "")){
          op_set bathy_set = cells;
          if (new_format) {
            bathy_nodes = op_decl_set_hdf5(filename_h5, "bathy_nodes");
            lifted_cells = op_decl_set_hdf5(filename_h5, "liftedCells");
            liftedcellsToBathyNodes = op_decl_map_hdf5(lifted_cells, bathy_nodes, N_NODESPERCELL,
                                                 filename_h5,
                                                 "liftedcellsToBathynodes");
            liftedcellsToCells = op_decl_map_hdf5(lifted_cells, cells, 1,
                                                 filename_h5,
                                                 "liftedcellsToCells");
            bathy_xy = op_decl_dat_hdf5(bathy_nodes, MESH_DIM, "float",
                                        filename_h5,
                                        "bathy_xy");
            initial_zb = op_decl_dat_hdf5(cells, 1, "float",
                                        filename_h5,
                                        "initial_zb");
            bathy_set = bathy_nodes;
          }

          // If one initBathymetry file is used
          if (strstr(events[i].streamName.c_str(), "%i") == NULL){
            n_initBathymetry = 1;
            temp_initBathymetry = (op_dat*) malloc(sizeof(op_dat));
            temp_initBathymetry[0] = op_decl_dat_hdf5(bathy_set, 1, "float",
                          filename_h5,
                          "initBathymetry");
          // If more initBathymetry files are used
          } else{
            if(timers[i].iend != INT_MAX) {
              n_initBathymetry = (timers[i].iend - timers[i].istart) / timers[i].istep + 1;
            } else {
              int tmp_iend = ftime/dtmax;
              n_initBathymetry = (tmp_iend-timers[i].istart)/timers[i].istep + 1;
            }
            op_printf("Reading %d consecutive InitBathymetry data arrays... ", n_initBathymetry);
            temp_initBathymetry = (op_dat*) malloc(n_initBathymetry * sizeof(op_dat));
            for(int k=0; k<n_initBathymetry; k++) {
                char dat_name[255];
                // iniBathymetry data is stored with sequential numbering instead of iteration step numbering!
                sprintf(dat_name,"initBathymetry%d",k);
                temp_initBathymetry[k] = op_decl_dat_hdf5(bathy_set, 1, "float",
                                filename_h5,
                                dat_name);
            }
            op_printf("done.\n");
          }
        }
      }
  }

  for (unsigned int i = 0; i < events.size(); i++) {
    if (!strcmp(events[i].className.c_str(), "InitBathyRelative")) {
      if (!new_format || n_initBathymetry < 1) {
        printf("Error, trying to use InitBathyRelative with new-format or without an initial InitBathymetry event\n");
        exit(-1);
      }
      initial_zb = temp_initBathymetry[0];
    }
  }


  if (op_is_root()) op_diagnostic_output();

  /*
   * Partitioning
   */
//  op_partition("PARMETIS", "GEOM", NULL, NULL, cellCenters);
//  op_partition("PTSCOTCH", "GEOM", NULL, NULL, cellCenters);
//  op_partition("", "", NULL, NULL, NULL);
  op_partition("PARMETIS", "KWAY", NULL, edgesToCells, NULL);
//  op_partition("PTSCOTCH", "KWAY", NULL, edgesToCells, NULL);
//  op_partition("PARMETIS", "GEOMKWAY", edges, edgesToCells, cellCenters);
//  op_partition("PARMETIS", "KWAY", NULL, NULL, NULL);
//  op_partition("PARMETIS", "KWAY", edges, edgesToCells, cellCenters);
//  op_partition("PTSCOTCH", "KWAY", NULL, cellsToEdges, NULL);
//  op_partition("PTSCOTCH", "KWAY", NULL, edgesToCells, NULL);
//  op_partition("PTSCOTCH", "KWAY", NULL, cellsToEdges, NULL);

  // Timer variables
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timers(&cpu_t1, &wall_t1);

  float *tmp_elem = NULL;
  if (num_outputLocation)
    outputLocation_dat = op_decl_dat_temp(outputLocation, 5, "float",
                                        tmp_elem,"outputLocation_dat");
  float zmin;
  op_dat z_zero = op_decl_dat_temp(cells, 1, "float",tmp_elem,"z_zero");
  //Very first Init loop
  processEvents(&timers, &events, 1/*firstTime*/, 1/*update timers*/, 0.0/*=dt*/, 1/*remove finished events*/, 2/*init loop, not pre/post*/, cells, values, cellVolumes, cellCenters, nodeCoords, cellsToNodes, temp_initEta, temp_initU, temp_initV, bathy_nodes, lifted_cells, liftedcellsToBathyNodes, liftedcellsToCells, bathy_xy, initial_zb, temp_initBathymetry, z_zero, n_initBathymetry, &zmin, outputLocation_map, outputLocation_dat, writeOption);
  /*
   *  Declaring temporary dats
  */
  op_dat values_new = op_decl_dat_temp(cells, 4, "float",tmp_elem,"values_new"); //tmp - cells - dim 4
  op_dat GradientatCell = op_decl_dat_temp(cells, 8, "float", tmp_elem, "GradientatCell");
  //SpaceDiscretization
  op_dat bathySource = op_decl_dat_temp(edges, 4, "float", tmp_elem, "bathySource"); //temp - edges - dim 2 (left & right)
  op_dat edgeFluxes = op_decl_dat_temp(edges, 3, "float", tmp_elem, "edgeFluxes"); //temp - edges - dim 4
  //NumericalFluxes
  op_dat maxEdgeEigenvalues = op_decl_dat_temp(edges, 1, "float", tmp_elem, "maxEdgeEigenvalues"); //temp - edges - dim 1
  //EvolveValuesRK22
  op_dat Lw_n = op_decl_dat_temp(cells, 4, "float", tmp_elem, "Lw_n"); //temp - cells - dim 4
  op_dat Lw_1 = op_decl_dat_temp(cells, 4, "float", tmp_elem, "Lw_1"); //temp - cells - dim 4
  op_dat w_1 = op_decl_dat_temp(cells, 4, "float", tmp_elem, "w_1"); //temp - cells - dim 4
  // q contains the max and min values of the physical variables surrounding each cell
  op_dat q = op_decl_dat_temp(cells, 8, "float", tmp_elem, "q"); //temp - cells - dim 8
  // lim is the limiter value for each physical variable defined on each cell
  op_dat lim = op_decl_dat_temp(cells, 4, "float", tmp_elem, "lim"); //temp - cells - dim 4
  double timestep;
  while (timestamp < ftime) {
		//process post_update==false events (usually Init events)
    processEvents(&timers, &events, 0, 0, 0.0, 0, 0, cells, values, cellVolumes, cellCenters, nodeCoords, cellsToNodes, temp_initEta, temp_initU, temp_initV, bathy_nodes,  lifted_cells, liftedcellsToBathyNodes, liftedcellsToCells, bathy_xy, initial_zb, temp_initBathymetry, z_zero, n_initBathymetry, &zmin, outputLocation_map, outputLocation_dat, writeOption);


#ifdef DEBUG
#endif
    {
      float minTimestep = 0.0;
      spaceDiscretization(values, Lw_n, &minTimestep,
          bathySource, edgeFluxes, maxEdgeEigenvalues,
          edgeNormals, edgeLength, cellVolumes, isBoundary,
          cells, edges, edgesToCells, cellsToEdges, cellsToCells, edgeCenters, cellCenters, GradientatCell, q, lim, &timestamp);
#ifdef DEBUG
      printf("Return of SpaceDiscretization #1 midPointConservative H %g U %g V %g Zb %g  \n", normcomp(Lw_n, 0), normcomp(Lw_n, 1),normcomp(Lw_n, 2),normcomp(Lw_n, 3));
#endif
      float dT = CFL * minTimestep;
      dT= dT < dtmax ? dT : dtmax;

      op_par_loop(EvolveValuesRK2_1, "EvolveValuesRK2_1", cells,
          op_arg_gbl(&dT,1,"float", OP_READ),
          op_arg_dat(Lw_n, -1, OP_ID, 4, "float", OP_READ),
          op_arg_dat(values, -1, OP_ID, 4, "float", OP_READ),
          op_arg_dat(w_1, -1, OP_ID, 4, "float", OP_WRITE));
#ifdef DEBUG
      printf("Return of SpaceDiscretization #1 midPointConservative H %g U %g V %g Zb %g  \n", normcomp(w_1, 0), normcomp(w_1, 1),normcomp(w_1, 2),normcomp(w_1, 3));
#endif

      float dummy = 0.0;
      spaceDiscretization(w_1, Lw_1, &dummy,
          bathySource, edgeFluxes, maxEdgeEigenvalues,
          edgeNormals, edgeLength, cellVolumes, isBoundary,
          cells, edges, edgesToCells, cellsToEdges,
          cellsToCells, edgeCenters, cellCenters, GradientatCell, q, lim, &timestamp);


      op_par_loop(EvolveValuesRK2_2, "EvolveValuesRK2_2", cells,
          op_arg_gbl(&dT,1,"float", OP_READ),
          op_arg_dat(Lw_1, -1, OP_ID, 4, "float", OP_READ),
          op_arg_dat(values, -1, OP_ID, 4, "float", OP_READ),
          op_arg_dat(w_1, -1, OP_ID, 4, "float", OP_READ),
          op_arg_dat(values_new, -1, OP_ID, 4, "float", OP_WRITE));


      timestep=dT;
      float Mn = 0.025f;
      op_par_loop(Friction_manning, "Friction_manning", cells,
          op_arg_gbl(&dT,1,"float", OP_READ),
          op_arg_gbl(&Mn,1,"float", OP_READ),
          op_arg_dat(values_new, -1, OP_ID, 4, "float", OP_RW));
    }
    op_par_loop(simulation_1, "simulation_1", cells,
        op_arg_dat(values, -1, OP_ID, 4, "float", OP_WRITE),
        op_arg_dat(values_new, -1, OP_ID, 4, "float", OP_READ));
//         printf("Return of SpaceDiscretization #1 midPointConservative H %g U %g V %g Zb %g  \n", normcomp(values_new, 0), normcomp(values_new, 1),normcomp(values_new, 2),normcomp(values_new, 3));

    itercount++;
    timestamp += timestep;
    processEvents(&timers, &events, 0, 1, timestep, 1, 1, cells, values, cellVolumes, cellCenters, nodeCoords, cellsToNodes,temp_initEta, temp_initU, temp_initV, bathy_nodes,    lifted_cells, liftedcellsToBathyNodes, liftedcellsToCells, bathy_xy, initial_zb, temp_initBathymetry, z_zero, n_initBathymetry, &zmin, outputLocation_map, outputLocation_dat, writeOption);
  }

  op_timers(&cpu_t2, &wall_t2);
  op_timing_output();
  op_printf("Main simulation runtime = \n%lf\n",wall_t2-wall_t1);


  if(op_is_root()) {
    int compressed = 0;
    if (locationData.n_points>0) {
      int len = locationData.time[0].size();
      int same = 1;
      for (int i = 1; i < locationData.n_points; i++) {
        if (locationData.time[0].size() != locationData.time[i].size()) same = 0;
      }
      if (same) compressed = 1;
    }
    if (compressed) {
      int len = locationData.time[0].size();
      int pts = locationData.n_points;
      int pts1 = locationData.n_points+1;
      //float *loc_data = (float*)malloc((locationData.n_points+1)*len*sizeof(float));
      float *loc_data = (float*)malloc((pts1*3)*len*sizeof(float));
      for (int i = 0; i < len; i++) {
        loc_data[i*pts1] = locationData.time[0][i];
        for (int j = 0; j < pts; j++) {
          loc_data[i*pts1+1+j] = locationData.value[j][i];
        }
      }
       /*loop for storing U*/
      for (int i = 0; i < len; i++) {
        loc_data[i*pts1+pts1*len] = locationData.time[0][i];
        for (int j = 0; j < pts; j++) {
          loc_data[i*pts1+pts1*len+1+j] = locationData.allvalues[j][4*i+1];
        }
      }
      	/*loop for storing V*/
      for (int i = 0; i < len; i++) {
        loc_data[i*pts1+2*pts1*len] = locationData.time[0][i];
        for (int j = 0; j < pts; j++) {
          loc_data[i*pts1+2*pts1*len+1+j] = locationData.allvalues[j][4*i+2];
        }
      }
      write_locations_hdf5(loc_data, pts1,len, "gauges.h5");
      write_locations_hdf5(loc_data+pts1*len, pts1,len, "gauges_U.h5");
      write_locations_hdf5(loc_data+2*pts1*len, pts1,len, "gauges_V.h5");
    } else {
      for (int i = 0; i < locationData.n_points; i++) {
        FILE* fp;
        fp = fopen(locationData.filename[i].c_str(), "w");
        if(fp == NULL) {
          op_printf("can't open file for write %s\n",locationData.filename[i].c_str());
          exit(-1);
        }
        for(unsigned int j=0; j<locationData.time[i].size() ; j++) {
          fprintf(fp, "%1.10f  %10.20g\n", locationData.time[i][j], locationData.value[i][j]);
        }
        if(fclose(fp)) {
          op_printf("can't close file %s\n",locationData.filename[i].c_str());
          exit(-1);
        }
      }
    }
    locationData.filename.clear();
    locationData.time.clear();
    locationData.value.clear();
    locationData.allvalues.clear();
  }

  for (int i = 0; i < timers.size(); i++) {
    if (timers[i].step == -1 && strcmp(events[i].className.c_str(), "OutputMaxElevation") == 0) {
      strcpy((char*)currentMaxElevation->name, "values");
      OutputSimulation(writeOption, &events[i], &timers[i], nodeCoords, cellsToNodes, currentMaxElevation, cells, &zmin);
      strcpy((char*)currentMaxElevation->name, "maxElevation");
    }
  }
  for (int i = 0; i < timers.size(); i++) {
    if (timers[i].step == -1 && strcmp(events[i].className.c_str(), "OutputMaxSpeed") == 0) {
      strcpy((char*)currentMaxSpeed->name, "values");
      OutputSimulation(writeOption, &events[i], &timers[i], nodeCoords, cellsToNodes, currentMaxSpeed, cells, &zmin);
      strcpy((char*)currentMaxSpeed->name, "maxSpeed");
    }
  }
  /*
   *	 Free temporary dats
   */
  //simulation
  if (op_free_dat_temp(values_new) < 0)
    op_printf("Error: temporary op_dat %s cannot be removed\n",values_new->name);
  if (op_free_dat_temp(GradientatCell) < 0)
    op_printf("Error: temporary op_dat %s cannot be removed\n",GradientatCell->name);

  //EvolveValuesRK2
  /*if (op_free_dat_temp(midPointConservative) < 0)
    op_printf("Error: temporary op_dat %s cannot be removed\n",midPointConservative->name);
  if (op_free_dat_temp(outConservative) < 0)
    op_printf("Error: temporary op_dat %s cannot be removed\n",outConservative->name);
  if (op_free_dat_temp(midPointConservative3) < 0)
    op_printf("Error: temporary op_dat %s cannot be removed\n",midPointConservative3->name);

  if (op_free_dat_temp(inConservative) < 0)
    op_printf("Error: temporary op_dat %s cannot be removed\n",inConservative->name);
  */
  if (op_free_dat_temp(Lw_n) < 0)
    op_printf("Error: temporary op_dat %s cannot be removed\n",Lw_n->name);
  if (op_free_dat_temp(Lw_1) < 0)
    op_printf("Error: temporary op_dat %s cannot be removed\n",Lw_1->name);
  if (op_free_dat_temp(w_1) < 0)
    op_printf("Error: temporary op_dat %s cannot be removed\n",w_1->name);
  //SpaceDiscretization
  if (op_free_dat_temp(bathySource) < 0)
    op_printf("Error: temporary op_dat %s cannot be removed\n",bathySource->name);
  if (op_free_dat_temp(edgeFluxes) < 0)
    op_printf("Error: temporary op_dat %s cannot be removed\n",edgeFluxes->name);
  //NumericalFluxes
  if (op_free_dat_temp(maxEdgeEigenvalues) < 0)
    op_printf("Error: temporary op_dat %s cannot be removed\n",maxEdgeEigenvalues->name);

  if (op_free_dat_temp(q)< 0)
    op_printf("Error: temporary op_dat %s cannot be removed\n",q->name);
  if (op_free_dat_temp(lim)< 0)
    op_printf("Error: temporary op_dat %s cannot be removed\n",lim->name);
  op_timers(&cpu_t2, &wall_t2);
  op_timing_output();
  op_printf("Max total runtime = \n%lf\n",wall_t2-wall_t1);

  op_exit();

  return 0;
}

#include "volna_common.h"
#include "op_lib_cpp.h"

void __check_hdf5_error(herr_t err, const char *file, const int line);

int timer_happens(TimerParams *p) {
  int result;
  result = ( p->t <= p->end && p->t >= p->start &&
      p->iter <= p->iend && p->iter >= p->istart );
  result =  result &&
      ( ( p->iter == 0) ||
          ( p->localIter == p->istep ) ||
          ( p->localTime >= p->step ) || (p->step == -1) );
  return result;
}

void write_locations_hdf5(float *data, int dim_cont, int dim_stride, const char *filename) {
  hid_t    file_id, dataset_id, dataspace_id; /* identifiers */
  hid_t    plist_id; 

  size_t   nelmts;
  unsigned flags, filter_info;
  H5Z_filter_t filter_type;

  herr_t   status;
  hsize_t  dims[1];
  hsize_t  cdims[1];

  int      idx;
  int      i,j, numfilt;

  /* Uncomment these variables to use SZIP compression 
  unsigned szip_options_mask;
  unsigned szip_pixels_per_block;
  */

  /* Create a file.  */
  file_id = H5Fcreate (filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  /* Write dimensions */
  dims[0] = 2;
  int dimensions[2] = {dim_cont, dim_stride};
  dataspace_id = H5Screate_simple(1, dims, NULL);
  dataset_id = H5Dcreate(file_id, "dims", H5T_NATIVE_INT, dataspace_id,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  //write data
  H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, dataspace_id, H5P_DEFAULT, (char*)dimensions);
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
  
  /* Create dataset "Compressed Data" in the group using absolute name.  */
  dims[0] = dim_cont*dim_stride;
  dataspace_id = H5Screate_simple (1, dims, NULL);

  plist_id  = H5Pcreate (H5P_DATASET_CREATE);

  /* Dataset must be chunked for compression */
  cdims[0] = dim_cont;
  status = H5Pset_chunk (plist_id, 1, cdims);

  /* Set ZLIB / DEFLATE Compression using compression level 6.
   * To use SZIP Compression comment out these lines. 
  */ 
  status = H5Pset_deflate (plist_id, 6); 

  /* Uncomment these lines to set SZIP Compression 
  szip_options_mask = H5_SZIP_NN_OPTION_MASK;
  szip_pixels_per_block = 16;
  status = H5Pset_szip (plist_id, szip_options_mask, szip_pixels_per_block);
  */
  
  dataset_id = H5Dcreate2 (file_id, "gauges", H5T_IEEE_F32BE, 
                          dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT); 


  status = H5Dwrite (dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

  status = H5Sclose (dataspace_id);
  status = H5Dclose (dataset_id);
  status = H5Pclose (plist_id);
  status = H5Fclose (file_id);
}

void read_events_hdf5(hid_t h5file, int num_events, std::vector<TimerParams> *timers, std::vector<EventParams> *events, int *num_outputLocation) {
  std::vector<float> timer_start(num_events);
  std::vector<float> timer_end(num_events);
  std::vector<float> timer_step(num_events);
  std::vector<int> timer_istart(num_events);
  std::vector<int> timer_iend(num_events);
  std::vector<int> timer_istep(num_events);

  std::vector<float> event_location_x(num_events);
  std::vector<float> event_location_y(num_events);
  std::vector<int> event_post_update(num_events);

  //const hsize_t num_events_hsize = num_events;
  check_hdf5_error(H5LTread_dataset(h5file, "timer_start", H5T_NATIVE_FLOAT, &timer_start[0]));
  check_hdf5_error(H5LTread_dataset(h5file, "timer_end", H5T_NATIVE_FLOAT, &timer_end[0]));
  check_hdf5_error(H5LTread_dataset(h5file, "timer_step", H5T_NATIVE_FLOAT, &timer_step[0]));
  check_hdf5_error(H5LTread_dataset(h5file, "timer_istart", H5T_NATIVE_INT, &timer_istart[0]));
  check_hdf5_error(H5LTread_dataset(h5file, "timer_iend", H5T_NATIVE_INT, &timer_iend[0]));
  check_hdf5_error(H5LTread_dataset(h5file, "timer_istep", H5T_NATIVE_INT, &timer_istep[0]));

  check_hdf5_error(H5LTread_dataset(h5file, "event_location_x", H5T_NATIVE_FLOAT, &event_location_x[0]));
  check_hdf5_error(H5LTread_dataset(h5file, "event_location_y", H5T_NATIVE_FLOAT, &event_location_y[0]));
  check_hdf5_error(H5LTread_dataset(h5file, "event_post_update", H5T_NATIVE_INT, &event_post_update[0]));

  /*
   * Convert Arrays to AoS
   */
  char buffer[22];
  std::vector<char> eventBuffer;
  int length = 0;
  for (int i = 0; i < num_events; i++) {
    (*timers)[i].start = timer_start[i];
    (*timers)[i].end = timer_end[i];
    (*timers)[i].step = timer_step[i];
    (*timers)[i].istart = timer_istart[i];
    (*timers)[i].iend = timer_iend[i];
    (*timers)[i].istep = timer_istep[i];

    (*events)[i].location_x = event_location_x[i];
    (*events)[i].location_y = event_location_y[i];
    (*events)[i].post_update = event_post_update[i];

    /*
     * If string can not handle a variable size char*, then use the commented lines
     */
    memset(buffer,0,22);
    sprintf(buffer, "event_className%d",i);
    check_hdf5_error(H5LTget_attribute_int(h5file, buffer, "length", &length));
    eventBuffer.resize(length);
    check_hdf5_error(H5LTread_dataset_string(h5file, buffer, &eventBuffer[0]));
    (*events)[i].className.assign(&eventBuffer[0], length);

    if (strcmp((*events)[i].className.c_str(), "OutputLocation") == 0)
      (*num_outputLocation)++;

    memset(buffer,0,22);
    sprintf(buffer, "event_formula%d",i);
    check_hdf5_error(H5LTget_attribute_int(h5file, buffer, "length", &length));
    eventBuffer.resize(length);
    check_hdf5_error(H5LTread_dataset_string(h5file, buffer, &eventBuffer[0]));
    (*events)[i].formula.assign(&eventBuffer[0], length);

    memset(buffer,0,22);
    sprintf(buffer, "event_streamName%d",i);
    check_hdf5_error(H5LTget_attribute_int(h5file, buffer, "length", &length));
    eventBuffer.resize(length);
    check_hdf5_error(H5LTread_dataset_string(h5file, buffer, &eventBuffer[0]));
    (*events)[i].streamName.assign(&eventBuffer[0], length);
  }
}

void processEvents(std::vector<TimerParams> *timers, std::vector<EventParams> *events, int firstTime, int updateTimers,
 									 float timeIncrement, int removeFinished, int initPrePost, op_set cells, op_dat values, op_dat cellVolumes,
									 op_dat cellCenters, op_dat nodeCoords, op_map cellsToNodes, op_dat temp_initEta, op_set bathy_nodes, op_map cellsToBathyNodes, op_dat bathy_xy, op_dat initial_zb, 
                   op_dat* temp_initBathymetry, int n_initBathymetry, BoreParams bore_params, GaussianLandslideParams gaussian_landslide_params, op_map outputLocation_map,
									 op_dat outputLocation_dat, int writeOption) {
  //  op_printf("processEvents()... \n");
  int size = (*timers).size();
  int i = 0;
  int j = 0;
  if (firstTime) {
    int k = 0;
    while (i < size) {
      if (strcmp((*events)[i].className.c_str(), "OutputLocation") == 0) {
        locationData.filename.push_back((*events)[i].streamName);
	(*events)[i].loc_index = k++;
      }
      i++;
    }
    locationData.n_points = locationData.filename.size();
    locationData.time.resize(locationData.n_points);
    locationData.value.resize(locationData.n_points);
    locationData.tmp = (float*) malloc(locationData.n_points*sizeof(float));
    i = 0;
  }
  while (i < size){
    if (timer_happens(&(*timers)[i]) && (initPrePost==2 || (*events)[i].post_update==initPrePost)) {
      if (strcmp((*events)[i].className.c_str(), "InitEta") == 0) {
        InitEta(cells, cellCenters, values, temp_initEta, temp_initEta!=NULL);
      } else if (strcmp((*events)[i].className.c_str(), "InitU") == 0) {
        InitU(cells, cellCenters, values);
      } else if (strcmp((*events)[i].className.c_str(), "InitV") == 0) {
        InitV(cells, cellCenters, values);
      } else if (strcmp((*events)[i].className.c_str(), "InitBathymetry") == 0) {
        // If initBathymetry is given by a formula (n_initBathymetry is 0), run InitBathymetry for formula
        if(n_initBathymetry == 0) {
          InitBathymetry(cells, cellCenters, values, NULL, 0, firstTime, bathy_nodes, cellsToBathyNodes, bathy_xy, initial_zb);
        }
        // If initBathymetry is given by 1 file, run InitBathymetry for that particular file
        if(n_initBathymetry == 1 ) {
          InitBathymetry(cells, cellCenters, values, *temp_initBathymetry, 1, firstTime, bathy_nodes, cellsToBathyNodes, bathy_xy, initial_zb);
        // Else if initBathymetry is given by multiple files, run InitBathymetry for those files
        } else if (n_initBathymetry > 1) {
          int k = ((*timers)[i].iter - (*timers)[i].istart) / (*timers)[i].istep;
          // Handle the case when InitBathymetry files are out for further bathymetry initalization: remove the event
          if(strcmp((*events)[i].className.c_str(), "InitBathymetry") == 0 && k<n_initBathymetry) {
            InitBathymetry(cells, cellCenters, values, temp_initBathymetry[k], 1, firstTime, bathy_nodes, cellsToBathyNodes, bathy_xy, initial_zb);
          }
        }
      } else if (strcmp((*events)[i].className.c_str(), "InitBore") == 0) {
        InitBore(cells, cellCenters, values, bore_params);
      } else if (strcmp((*events)[i].className.c_str(), "InitGaussianLandslide") == 0) {
        InitGaussianLandslide(cells, cellCenters, values, gaussian_landslide_params, firstTime);
      } else if (strcmp((*events)[i].className.c_str(), "OutputTime") == 0) {
        OutputTime(&(*timers)[i]);
        //op_printf("Output iter: %d \n", (*timers)[i].iter);
      } else if (strcmp((*events)[i].className.c_str(), "OutputConservedQuantities") == 0) {
        OutputConservedQuantities(cells, cellVolumes, values);
      } else if (strcmp((*events)[i].className.c_str(), "OutputLocation") == 0) {
        OutputLocation(&(*events)[i], j, &(*timers)[i], cells, nodeCoords, cellsToNodes, values, outputLocation_map, outputLocation_dat);
				    j++;
      } else if (strcmp((*events)[i].className.c_str(), "OutputSimulation") == 0) {
        // Remove comment if needed:
        // 0 - HDF5 output
        // 1 - VTK ASCII output
        // 2 - VTK Binary output
        OutputSimulation(writeOption, &(*events)[i], &(*timers)[i], nodeCoords, cellsToNodes, values);
      } else if (strcmp((*events)[i].className.c_str(), "OutputMaxElevation") == 0) {
        // 0 - HDF5 output
        // 1 - VTK ASCII output
        // 2 - VTK Binary output
        OutputMaxElevation(writeOption, &(*events)[i], &(*timers)[i], nodeCoords, cellsToNodes, values, cells);
      } else {
        op_printf("Unrecognized event %s\n", (*events)[i].className.c_str());
        exit(-1);
      }
      //timer.LocalReset();
      (*timers)[i].localIter = 0;
      (*timers)[i].localTime = 0;
    }
    if (updateTimers) {
      //op_printf("if(updatesTimers) => true \n");
      //timer.update()
      (*timers)[i].t+= timeIncrement;
      (*timers)[i].iter += 1;
      (*timers)[i].localIter += 1;
      (*timers)[i].localTime += timeIncrement;
    }

    if (removeFinished) {
      //op_printf("if(removeFinished) => true \n");
      //Remove finished events
      if (((*timers)[i].iter >= (*timers)[i].iend) || ((*timers)[i].t >= (*timers)[i].end)) {
        (*timers).erase((*timers).begin()+i);
        (*events).erase((*events).begin()+i);
        size--;
      } else i++;
    } else i++;
  }
}

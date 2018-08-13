//
// auto-generated by op2.py
//

//user function
inline void initBathymetry_large(float *values, const float *cellCenter,
 const float *node0, const float *node1, const float *node2,
 const float *bathy0, const float *bathy1, const float *bathy2) {


 bool isInside = false;



 float xmin = MIN(MIN(node0[0], node1[0]), node2[0]);
 float xmax = MAX(MAX(node0[0], node1[0]), node2[0]);
 float ymin = MIN(MIN(node0[1], node1[1]), node2[1]);
 float ymax = MAX(MAX(node0[1], node1[1]), node2[1]);

 if ( ( cellCenter[0] < xmin ) || ( cellCenter[0] > xmax ) ||
  ( cellCenter[1] < ymin ) || ( cellCenter[1] > ymax ) ) {
  isInside = false;
 }else{



  float insider = 1.0f;
  float p[2] = {cellCenter[0], cellCenter[1]};


    if ( (node0[0] - node2[0]) * (node1[1] - node2[1]) - (node0[1] - node2[1]) * (node1[0] - node2[0]) > 0 ) {
     insider = (node0[0] - node2[0]) * (p[1] - node2[1]) - (node0[1] - node2[1]) * (p[0] - node2[0]);
     insider *= (node0[0] - p[0]) * (node1[1] - p[1]) - (node0[1] - p[1]) * (node1[0] - p[0]);
     insider *= (node1[0] - p[0]) * (node2[1] - p[1]) - (node1[1] - p[1]) * (node2[0] - p[0]);
    }
    else {
     insider = (node0[0] - node1[0]) * (p[1] - node1[1]) - (node0[1] - node1[1]) * (p[0] - node1[0]);
     insider *= (node0[0] - p[0]) * (node2[1] - p[1]) - (node0[1] - p[1]) * (node2[0] - p[0]);
     insider *= (node2[0] - p[0]) * (node1[1] - p[1]) - (node2[1] - p[1]) * (node1[0] - p[0]);
    }
    isInside = insider >= 0.0f;
 }


  if (isInside) {

    float a = (node1[1]-node0[1])*(*bathy2-*bathy0)-(node2[1]-node0[1])*(*bathy1-*bathy0);
    float b = -(node1[0]-node0[0])*(*bathy2-*bathy0)+(node2[0]-node0[0])*(*bathy1-*bathy0);
    float c = (node1[0]-node0[0])*(node2[1]-node0[1])-(node2[0]-node0[0])*(node1[1]-node0[1]);

    values[3] += *bathy0 - (a*(cellCenter[0]-node0[0]) + b*(cellCenter[1]-node0[1]))/c;
  }
}
#ifdef VECTORIZE
//user function -- modified for vectorisation
void initBathymetry_large_vec( float values[*][SIMD_VEC], const float cellCenter[*][SIMD_VEC], const float node0[*][SIMD_VEC], const float node1[*][SIMD_VEC], const float node2[*][SIMD_VEC], const float bathy0[*][SIMD_VEC], const float bathy1[*][SIMD_VEC], const float bathy2[*][SIMD_VEC], int idx ) {

 bool isInside = false;


 float xmin = MIN(MIN(node0[0][idx], node1[0][idx]), node2[0][idx]);
 float xmax = MAX(MAX(node0[0][idx], node1[0][idx]), node2[0][idx]);
 float ymin = MIN(MIN(node0[1][idx], node1[1][idx]), node2[1][idx]);
 float ymax = MAX(MAX(node0[1][idx], node1[1][idx]), node2[1][idx]);

 if ( ( cellCenter[0][idx] < xmin ) || ( cellCenter[0][idx] > xmax ) ||
  ( cellCenter[1][idx] < ymin ) || ( cellCenter[1][idx] > ymax ) ) {
  isInside = false;
 }else{


  float insider = 1.0f;
  float p[2] = {cellCenter[0][idx], cellCenter[1][idx]};

    if ( (node0[0][idx] - node2[0][idx]) * (node1[1][idx] - node2[1][idx]) - (node0[1][idx] - node2[1][idx]) * (node1[0][idx] - node2[0][idx]) > 0 ) {
     insider = (node0[0][idx] - node2[0][idx]) * (p[1] - node2[1][idx]) - (node0[1][idx] - node2[1][idx]) * (p[0] - node2[0][idx]);
     insider *= (node0[0][idx] - p[0]) * (node1[1][idx] - p[1]) - (node0[1][idx] - p[1]) * (node1[0][idx] - p[0]);
     insider *= (node1[0][idx] - p[0]) * (node2[1][idx] - p[1]) - (node1[1][idx] - p[1]) * (node2[0][idx] - p[0]);
    }
    else {
     insider = (node0[0][idx] - node1[0][idx]) * (p[1] - node1[1][idx]) - (node0[1][idx] - node1[1][idx]) * (p[0] - node1[0][idx]);
     insider *= (node0[0][idx] - p[0]) * (node2[1][idx] - p[1]) - (node0[1][idx] - p[1]) * (node2[0][idx] - p[0]);
     insider *= (node2[0][idx] - p[0]) * (node1[1][idx] - p[1]) - (node2[1][idx] - p[1]) * (node1[0][idx] - p[0]);
    }
    isInside = insider >= 0.0f;
 }

  if (isInside) {

    float a = (node1[1][idx]-node0[1][idx])*(bathy2[0][idx]-bathy0[0][idx])-(node2[1][idx]-node0[1][idx])*(bathy1[0][idx]-bathy0[0][idx]);
    float b = -(node1[0][idx]-node0[0][idx])*(bathy2[0][idx]-bathy0[0][idx])+(node2[0][idx]-node0[0][idx])*(bathy1[0][idx]-bathy0[0][idx]);
    float c = (node1[0][idx]-node0[0][idx])*(node2[1][idx]-node0[1][idx])-(node2[0][idx]-node0[0][idx])*(node1[1][idx]-node0[1][idx]);

    values[3][idx] += bathy0[0][idx]- (a*(cellCenter[0][idx]-node0[0][idx]) + b*(cellCenter[1][idx]-node0[1][idx]))/c;
  }
}
#endif

// host stub function
void op_par_loop_initBathymetry_large(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7){

  int nargs = 8;
  op_arg args[8];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;
  args[7] = arg7;
  //create aligned pointers for dats
  ALIGNED_float       float * __restrict__ ptr0 = (float *) arg0.data;
  __assume_aligned(ptr0,float_ALIGN);
  ALIGNED_float const float * __restrict__ ptr1 = (float *) arg1.data;
  __assume_aligned(ptr1,float_ALIGN);
  ALIGNED_float const float * __restrict__ ptr2 = (float *) arg2.data;
  __assume_aligned(ptr2,float_ALIGN);
  ALIGNED_float const float * __restrict__ ptr3 = (float *) arg3.data;
  __assume_aligned(ptr3,float_ALIGN);
  ALIGNED_float const float * __restrict__ ptr4 = (float *) arg4.data;
  __assume_aligned(ptr4,float_ALIGN);
  ALIGNED_float const float * __restrict__ ptr5 = (float *) arg5.data;
  __assume_aligned(ptr5,float_ALIGN);
  ALIGNED_float const float * __restrict__ ptr6 = (float *) arg6.data;
  __assume_aligned(ptr6,float_ALIGN);
  ALIGNED_float const float * __restrict__ ptr7 = (float *) arg7.data;
  __assume_aligned(ptr7,float_ALIGN);

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(20);
  op_timers_core(&cpu_t1, &wall_t1);

  if (OP_diags>2) {
    printf(" kernel routine with indirection: initBathymetry_large\n");
  }

  int exec_size = op_mpi_halo_exchanges(set, nargs, args);

  if (exec_size >0) {

    #ifdef VECTORIZE
    #pragma novector
    for ( int n=0; n<(exec_size/SIMD_VEC)*SIMD_VEC; n+=SIMD_VEC ){
      if (n+SIMD_VEC >= set->core_size) {
        op_mpi_wait_all(nargs, args);
      }
      ALIGNED_float float dat0[4][SIMD_VEC];
      ALIGNED_float float dat1[2][SIMD_VEC];
      ALIGNED_float float dat2[2][SIMD_VEC];
      ALIGNED_float float dat3[2][SIMD_VEC];
      ALIGNED_float float dat4[2][SIMD_VEC];
      ALIGNED_float float dat5[1][SIMD_VEC];
      ALIGNED_float float dat6[1][SIMD_VEC];
      ALIGNED_float float dat7[1][SIMD_VEC];
      #pragma simd
      for ( int i=0; i<SIMD_VEC; i++ ){
        int idx1_2 = 2 * arg0.map_data[(n+i) * arg0.map->dim + 0];
        int idx2_2 = 2 * arg2.map_data[(n+i) * arg2.map->dim + 0];
        int idx3_2 = 2 * arg2.map_data[(n+i) * arg2.map->dim + 1];
        int idx4_2 = 2 * arg2.map_data[(n+i) * arg2.map->dim + 2];
        int idx5_1 = 1 * arg2.map_data[(n+i) * arg2.map->dim + 0];
        int idx6_1 = 1 * arg2.map_data[(n+i) * arg2.map->dim + 1];
        int idx7_1 = 1 * arg2.map_data[(n+i) * arg2.map->dim + 2];

        dat0[0][i] = 0.0;
        dat0[1][i] = 0.0;
        dat0[2][i] = 0.0;
        dat0[3][i] = 0.0;

        dat1[0][i] = (ptr1)[idx1_2 + 0];
        dat1[1][i] = (ptr1)[idx1_2 + 1];

        dat2[0][i] = (ptr2)[idx2_2 + 0];
        dat2[1][i] = (ptr2)[idx2_2 + 1];

        dat3[0][i] = (ptr3)[idx3_2 + 0];
        dat3[1][i] = (ptr3)[idx3_2 + 1];

        dat4[0][i] = (ptr4)[idx4_2 + 0];
        dat4[1][i] = (ptr4)[idx4_2 + 1];

        dat5[0][i] = (ptr5)[idx5_1 + 0];

        dat6[0][i] = (ptr6)[idx6_1 + 0];

        dat7[0][i] = (ptr7)[idx7_1 + 0];

      }
      #pragma simd
      for ( int i=0; i<SIMD_VEC; i++ ){
        initBathymetry_large_vec(
          dat0,
          dat1,
          dat2,
          dat3,
          dat4,
          dat5,
          dat6,
          dat7,
          i);
      }
      for ( int i=0; i<SIMD_VEC; i++ ){
        int idx0_4 = 4 * arg0.map_data[(n+i) * arg0.map->dim + 0];

        (ptr0)[idx0_4 + 0] += dat0[0][i];
        (ptr0)[idx0_4 + 1] += dat0[1][i];
        (ptr0)[idx0_4 + 2] += dat0[2][i];
        (ptr0)[idx0_4 + 3] += dat0[3][i];

      }
    }

    //remainder
    for ( int n=(exec_size/SIMD_VEC)*SIMD_VEC; n<exec_size; n++ ){
    #else
    for ( int n=0; n<exec_size; n++ ){
    #endif
      if (n==set->core_size) {
        op_mpi_wait_all(nargs, args);
      }
      int map0idx = arg0.map_data[n * arg0.map->dim + 0];
      int map2idx = arg2.map_data[n * arg2.map->dim + 0];
      int map3idx = arg2.map_data[n * arg2.map->dim + 1];
      int map4idx = arg2.map_data[n * arg2.map->dim + 2];

      initBathymetry_large(
        &(ptr0)[4 * map0idx],
        &(ptr1)[2 * map0idx],
        &(ptr2)[2 * map2idx],
        &(ptr3)[2 * map3idx],
        &(ptr4)[2 * map4idx],
        &(ptr5)[1 * map2idx],
        &(ptr6)[1 * map3idx],
        &(ptr7)[1 * map4idx]);
    }
  }

  if (exec_size == 0 || exec_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[20].name      = name;
  OP_kernels[20].count    += 1;
  OP_kernels[20].time     += wall_t2 - wall_t1;
  OP_kernels[20].transfer += (float)set->size * arg0.size * 2.0f;
  OP_kernels[20].transfer += (float)set->size * arg1.size;
  OP_kernels[20].transfer += (float)set->size * arg2.size;
  OP_kernels[20].transfer += (float)set->size * arg5.size;
  OP_kernels[20].transfer += (float)set->size * arg0.map->dim * 4.0f;
  OP_kernels[20].transfer += (float)set->size * arg2.map->dim * 4.0f;
}

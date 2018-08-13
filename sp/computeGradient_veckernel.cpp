//
// auto-generated by op2.py
//

//user function
inline void computeGradient(const float *center,
                            const float *neighbour1,
                            const float *neighbour2,
                            const float *neighbour3,
                            const float *cellCenter,
                            const float *nb1Center,
                            const float *nb2Center,
                            const float *nb3Center,
                            float *q, float *out) //OP_WRITE
{
  // Only reconstruct if the cell is not a touching the edge
  // Least-Squares Gradient Reconstruction
  if( cellCenter[0] != nb3Center[0] && cellCenter[1] != nb3Center[1]){
    float total, Rhs[8];
    float dh[3], dz[3],du[3], dv[3], weights[3];
    float Gram[2][2], inverse[2][2], delta[3][2];
    float x = cellCenter[0];
    float y = cellCenter[1];
    // Finding delta_x,delta_y for the neighbouring cells.
    delta[0][0] =  (nb1Center[0] - x);
    delta[0][1] =  (nb1Center[1] - y);

    delta[1][0] =  (nb2Center[0] - x);
    delta[1][1] =  (nb2Center[1] - y);

    delta[2][0] =  (nb3Center[0] - x);
    delta[2][1] =  (nb3Center[1] - y);
    // Calculating the weights coefficients based on the distance between
    // neighbouring cells and center cell.
    weights[0] = sqrt(delta[0][0] * delta[0][0] + delta[0][1] * delta[0][1]);
    weights[1] = sqrt(delta[1][0] * delta[1][0] + delta[1][1] * delta[1][1]);
    weights[2] = sqrt(delta[2][0] * delta[2][0] + delta[2][1] * delta[2][1]);

    total = weights[0] + weights[1] + weights[2];
    weights[0] = total/weights[0];
    weights[1] = total/weights[1];
    weights[2] = total/ weights[2];

    delta[0][0] *= weights[0];
    delta[0][1] *= weights[0];

    delta[1][0] *= weights[1];
    delta[1][1] *= weights[1];

    delta[2][0] *= weights[2];
    delta[2][1] *= weights[2];
    // Setting up the Gram matrix
    Gram[0][0] = ((delta[0][0]*delta[0][0]) + (delta[1][0] *delta[1][0]) + (delta[2][0] *delta[2][0]));
    Gram[0][1] = ((delta[0][0]*delta[0][1]) + (delta[1][0] *delta[1][1]) + (delta[2][0] *delta[2][1]));
    Gram[1][0] = ((delta[0][0]*delta[0][1]) + (delta[1][0] *delta[1][1]) + (delta[2][0] *delta[2][1]));
    Gram[1][1] = ((delta[0][1]*delta[0][1]) + (delta[1][1] *delta[1][1]) + (delta[2][1] *delta[2][1]));
    // Finding the inverse of the determinant
    float det = 1.0 / (Gram[0][0]*Gram[1][1] - Gram[0][1]*Gram[1][0]);
    inverse[0][0] = det * Gram[1][1];
    inverse[0][1] = det * (- Gram[0][1]);
    inverse[1][0] = det * (-Gram[1][0]);
    inverse[1][1] = det * Gram[0][0];
    // Setting up the RHS
    dh[0] = neighbour1[0] - center[0];
    dh[1] = neighbour2[0] - center[0];
    dh[2] = neighbour3[0] - center[0];
    dh[0] *= weights[0];
    dh[1] *= weights[1];
    dh[2] *= weights[2];

    dz[0] = neighbour1[3] - center[3];
    dz[1] = neighbour2[3] - center[3];
    dz[2] = neighbour3[3] - center[3];
    dz[0] *= weights[0];
    dz[1] *= weights[1];
    dz[2] *= weights[2];

    du[0] = neighbour1[1] - center[1];
    du[1] = neighbour2[1] - center[1];
    du[2] = neighbour3[1] - center[1];
    du[0] *= weights[0];
    du[1] *= weights[1];
    du[2] *= weights[2];

    dv[0] = neighbour1[2] - center[2];
    dv[1] = neighbour2[2] - center[2];
    dv[2] = neighbour3[2] - center[2];
    dv[0] *= weights[0];
    dv[1] *= weights[1];
    dv[2] *= weights[2];

    Rhs[0] = (delta[0][0]*dh[0]) + (delta[1][0]*dh[1]) + (delta[2][0]*dh[2]);
    Rhs[1] = (delta[0][1]*dh[0]) + (delta[1][1]*dh[1]) + (delta[2][1]*dh[2]);
    out[0] = (inverse[0][0] * Rhs[0]) + (inverse[0][1] * Rhs[1]);
    out[1] = (inverse[1][0] * Rhs[0]) + (inverse[1][1] * Rhs[1]);

    Rhs[2] = (delta[0][0]*du[0]) + (delta[1][0]*du[1]) + (delta[2][0]*du[2]);
    Rhs[3] = (delta[0][1]*du[0]) + (delta[1][1]*du[1]) + (delta[2][1]*du[2]);
    out[2] = (inverse[0][0] * Rhs[2]) + (inverse[0][1] * Rhs[3]);
    out[3] = (inverse[1][0] * Rhs[2]) + (inverse[1][1] * Rhs[3]);

    Rhs[4] = (delta[0][0]*dv[0]) + (delta[1][0]*dv[1]) + (delta[2][0]*dv[2]);
    Rhs[5] = (delta[0][1]*dv[0]) + (delta[1][1]*dv[1]) + (delta[2][1]*dv[2]);
    out[4] = (inverse[0][0] * Rhs[4]) + (inverse[0][1] * Rhs[5]);
    out[5] = (inverse[1][0] * Rhs[4]) + (inverse[1][1] * Rhs[5]);

    Rhs[6] = (delta[0][0]*dz[0]) + (delta[1][0]*dz[1]) + (delta[2][0]*dz[2]);
    Rhs[7] = (delta[0][1]*dz[0]) + (delta[1][1]*dz[1]) + (delta[2][1]*dz[2]);
    out[6] = (inverse[0][0] * Rhs[6]) + (inverse[0][1] * Rhs[7]);
    out[7] = (inverse[1][0] * Rhs[6]) + (inverse[1][1] * Rhs[7]);
 }else {
    // Gradients for the edge cells are set to zero.
    out[0] = 0.0f;
    out[1] = 0.0f;
    out[2] = 0.0f;
    out[3] = 0.0f;
    out[4] = 0.0f;
    out[5] = 0.0f;
    out[6] = 0.0f;
    out[7] = 0.0f;
 }
  // Computed the local max and min values for H,U,V,Z
  // q[0] - Hmin , q[1] - Hmax
  // q[2] - Umin , q[3] - Umax
  // q[4] - Vmin , q[5] - Vmax
  // q[6] - Zmin , q[7] - Zmax
  q[0] = center[0] < neighbour1[0] ? center[0] : neighbour1[0];
  q[0] = q[0] < neighbour2[0] ? q[0] : neighbour2[0];
  q[0] = q[0] < neighbour3[0] ? q[0] : neighbour3[0];
  q[1] = center[0] > neighbour1[0] ? center[0] : neighbour1[0];
  q[1] = q[1] > neighbour2[0] ? q[1] : neighbour2[0];
  q[1] = q[1] > neighbour3[0] ? q[1] : neighbour3[0];

  q[2] = center[1] < neighbour1[1] ? center[1] : neighbour1[1];
  q[2] = q[2] < neighbour2[1] ? q[2] : neighbour2[1];
  q[2] = q[2] < neighbour3[1] ? q[2] : neighbour3[1];
  q[3] = center[1] > neighbour1[1] ? center[1] : neighbour1[1];
  q[3] = q[3] > neighbour2[1] ? q[3] : neighbour2[1];
  q[3] = q[3] > neighbour3[1] ? q[3] : neighbour3[1];

  q[4] = center[2] < neighbour1[2] ? center[2] : neighbour1[2];
  q[4] = q[4] < neighbour2[2] ? q[4] : neighbour2[2];
  q[4] = q[4] < neighbour3[2] ? q[4] : neighbour3[2];
  q[5] = center[2] > neighbour1[2] ? center[2] : neighbour1[2];
  q[5] = q[5] > neighbour2[2] ? q[5] : neighbour2[2];
  q[5] = q[5] > neighbour3[2] ? q[5] : neighbour3[2];

  q[6] = center[3] < neighbour1[3] ? center[3] : neighbour1[3];
  q[6] = q[6] < neighbour2[3] ? q[6] : neighbour2[3];
  q[6] = q[6] < neighbour3[3] ? q[6] : neighbour3[3];
  q[7] = center[3] > neighbour1[3] ? center[3] : neighbour1[3];
  q[7] = q[7] > neighbour2[3] ? q[7] : neighbour2[3];
  q[7] = q[7] > neighbour3[3] ? q[7] : neighbour3[3];
}
#ifdef VECTORIZE
//user function -- modified for vectorisation
void computeGradient_vec( const float *center, const float neighbour1[*][SIMD_VEC], const float neighbour2[*][SIMD_VEC], const float neighbour3[*][SIMD_VEC], const float *cellCenter, const float nb1Center[*][SIMD_VEC], const float nb2Center[*][SIMD_VEC], const float nb3Center[*][SIMD_VEC], float *q, float *out, int idx ) {


  if( cellCenter[0] != nb3Center[0][idx] && cellCenter[1] != nb3Center[1][idx]){
    float total, Rhs[8];
    float dh[3], dz[3],du[3], dv[3], weights[3];
    float Gram[2][2], inverse[2][2], delta[3][2];
    float x = cellCenter[0];
    float y = cellCenter[1];

    delta[0][0] =  (nb1Center[0][idx] - x);
    delta[0][1] =  (nb1Center[1][idx] - y);

    delta[1][0] =  (nb2Center[0][idx] - x);
    delta[1][1] =  (nb2Center[1][idx] - y);

    delta[2][0] =  (nb3Center[0][idx] - x);
    delta[2][1] =  (nb3Center[1][idx] - y);


    weights[0] = sqrt(delta[0][0] * delta[0][0] + delta[0][1] * delta[0][1]);
    weights[1] = sqrt(delta[1][0] * delta[1][0] + delta[1][1] * delta[1][1]);
    weights[2] = sqrt(delta[2][0] * delta[2][0] + delta[2][1] * delta[2][1]);

    total = weights[0] + weights[1] + weights[2];
    weights[0] = total/weights[0];
    weights[1] = total/weights[1];
    weights[2] = total/ weights[2];

    delta[0][0] *= weights[0];
    delta[0][1] *= weights[0];

    delta[1][0] *= weights[1];
    delta[1][1] *= weights[1];

    delta[2][0] *= weights[2];
    delta[2][1] *= weights[2];

    Gram[0][0] = ((delta[0][0]*delta[0][0]) + (delta[1][0] *delta[1][0]) + (delta[2][0] *delta[2][0]));
    Gram[0][1] = ((delta[0][0]*delta[0][1]) + (delta[1][0] *delta[1][1]) + (delta[2][0] *delta[2][1]));
    Gram[1][0] = ((delta[0][0]*delta[0][1]) + (delta[1][0] *delta[1][1]) + (delta[2][0] *delta[2][1]));
    Gram[1][1] = ((delta[0][1]*delta[0][1]) + (delta[1][1] *delta[1][1]) + (delta[2][1] *delta[2][1]));

    float det = 1.0 / (Gram[0][0]*Gram[1][1] - Gram[0][1]*Gram[1][0]);
    inverse[0][0] = det * Gram[1][1];
    inverse[0][1] = det * (- Gram[0][1]);
    inverse[1][0] = det * (-Gram[1][0]);
    inverse[1][1] = det * Gram[0][0];

    dh[0] = neighbour1[0][idx] - center[0];
    dh[1] = neighbour2[0][idx] - center[0];
    dh[2] = neighbour3[0][idx] - center[0];
    dh[0] *= weights[0];
    dh[1] *= weights[1];
    dh[2] *= weights[2];

    dz[0] = neighbour1[3][idx] - center[3];
    dz[1] = neighbour2[3][idx] - center[3];
    dz[2] = neighbour3[3][idx] - center[3];
    dz[0] *= weights[0];
    dz[1] *= weights[1];
    dz[2] *= weights[2];

    du[0] = neighbour1[1][idx] - center[1];
    du[1] = neighbour2[1][idx] - center[1];
    du[2] = neighbour3[1][idx] - center[1];
    du[0] *= weights[0];
    du[1] *= weights[1];
    du[2] *= weights[2];

    dv[0] = neighbour1[2][idx] - center[2];
    dv[1] = neighbour2[2][idx] - center[2];
    dv[2] = neighbour3[2][idx] - center[2];
    dv[0] *= weights[0];
    dv[1] *= weights[1];
    dv[2] *= weights[2];

    Rhs[0] = (delta[0][0]*dh[0]) + (delta[1][0]*dh[1]) + (delta[2][0]*dh[2]);
    Rhs[1] = (delta[0][1]*dh[0]) + (delta[1][1]*dh[1]) + (delta[2][1]*dh[2]);
    out[0] = (inverse[0][0] * Rhs[0]) + (inverse[0][1] * Rhs[1]);
    out[1] = (inverse[1][0] * Rhs[0]) + (inverse[1][1] * Rhs[1]);

    Rhs[2] = (delta[0][0]*du[0]) + (delta[1][0]*du[1]) + (delta[2][0]*du[2]);
    Rhs[3] = (delta[0][1]*du[0]) + (delta[1][1]*du[1]) + (delta[2][1]*du[2]);
    out[2] = (inverse[0][0] * Rhs[2]) + (inverse[0][1] * Rhs[3]);
    out[3] = (inverse[1][0] * Rhs[2]) + (inverse[1][1] * Rhs[3]);

    Rhs[4] = (delta[0][0]*dv[0]) + (delta[1][0]*dv[1]) + (delta[2][0]*dv[2]);
    Rhs[5] = (delta[0][1]*dv[0]) + (delta[1][1]*dv[1]) + (delta[2][1]*dv[2]);
    out[4] = (inverse[0][0] * Rhs[4]) + (inverse[0][1] * Rhs[5]);
    out[5] = (inverse[1][0] * Rhs[4]) + (inverse[1][1] * Rhs[5]);

    Rhs[6] = (delta[0][0]*dz[0]) + (delta[1][0]*dz[1]) + (delta[2][0]*dz[2]);
    Rhs[7] = (delta[0][1]*dz[0]) + (delta[1][1]*dz[1]) + (delta[2][1]*dz[2]);
    out[6] = (inverse[0][0] * Rhs[6]) + (inverse[0][1] * Rhs[7]);
    out[7] = (inverse[1][0] * Rhs[6]) + (inverse[1][1] * Rhs[7]);
 }else {

    out[0] = 0.0f;
    out[1] = 0.0f;
    out[2] = 0.0f;
    out[3] = 0.0f;
    out[4] = 0.0f;
    out[5] = 0.0f;
    out[6] = 0.0f;
    out[7] = 0.0f;
 }





  q[0] = center[0] < neighbour1[0][idx] ? center[0] : neighbour1[0][idx];
  q[0] = q[0] < neighbour2[0][idx] ? q[0] : neighbour2[0][idx];
  q[0] = q[0] < neighbour3[0][idx] ? q[0] : neighbour3[0][idx];
  q[1] = center[0] > neighbour1[0][idx] ? center[0] : neighbour1[0][idx];
  q[1] = q[1] > neighbour2[0][idx] ? q[1] : neighbour2[0][idx];
  q[1] = q[1] > neighbour3[0][idx] ? q[1] : neighbour3[0][idx];

  q[2] = center[1] < neighbour1[1][idx] ? center[1] : neighbour1[1][idx];
  q[2] = q[2] < neighbour2[1][idx] ? q[2] : neighbour2[1][idx];
  q[2] = q[2] < neighbour3[1][idx] ? q[2] : neighbour3[1][idx];
  q[3] = center[1] > neighbour1[1][idx] ? center[1] : neighbour1[1][idx];
  q[3] = q[3] > neighbour2[1][idx] ? q[3] : neighbour2[1][idx];
  q[3] = q[3] > neighbour3[1][idx] ? q[3] : neighbour3[1][idx];

  q[4] = center[2] < neighbour1[2][idx] ? center[2] : neighbour1[2][idx];
  q[4] = q[4] < neighbour2[2][idx] ? q[4] : neighbour2[2][idx];
  q[4] = q[4] < neighbour3[2][idx] ? q[4] : neighbour3[2][idx];
  q[5] = center[2] > neighbour1[2][idx] ? center[2] : neighbour1[2][idx];
  q[5] = q[5] > neighbour2[2][idx] ? q[5] : neighbour2[2][idx];
  q[5] = q[5] > neighbour3[2][idx] ? q[5] : neighbour3[2][idx];

  q[6] = center[3] < neighbour1[3][idx] ? center[3] : neighbour1[3][idx];
  q[6] = q[6] < neighbour2[3][idx] ? q[6] : neighbour2[3][idx];
  q[6] = q[6] < neighbour3[3][idx] ? q[6] : neighbour3[3][idx];
  q[7] = center[3] > neighbour1[3][idx] ? center[3] : neighbour1[3][idx];
  q[7] = q[7] > neighbour2[3][idx] ? q[7] : neighbour2[3][idx];
  q[7] = q[7] > neighbour3[3][idx] ? q[7] : neighbour3[3][idx];
}
#endif

// host stub function
void op_par_loop_computeGradient(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9){

  int nargs = 10;
  op_arg args[10];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;
  args[7] = arg7;
  args[8] = arg8;
  args[9] = arg9;
  //create aligned pointers for dats
  ALIGNED_float const float * __restrict__ ptr0 = (float *) arg0.data;
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
  ALIGNED_float       float * __restrict__ ptr8 = (float *) arg8.data;
  __assume_aligned(ptr8,float_ALIGN);
  ALIGNED_float       float * __restrict__ ptr9 = (float *) arg9.data;
  __assume_aligned(ptr9,float_ALIGN);

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(5);
  op_timers_core(&cpu_t1, &wall_t1);

  if (OP_diags>2) {
    printf(" kernel routine with indirection: computeGradient\n");
  }

  int exec_size = op_mpi_halo_exchanges(set, nargs, args);

  if (exec_size >0) {

    #ifdef VECTORIZE
    #pragma novector
    for ( int n=0; n<(exec_size/SIMD_VEC)*SIMD_VEC; n+=SIMD_VEC ){
      if (n+SIMD_VEC >= set->core_size) {
        op_mpi_wait_all(nargs, args);
      }
      ALIGNED_float float dat1[4][SIMD_VEC];
      ALIGNED_float float dat2[4][SIMD_VEC];
      ALIGNED_float float dat3[4][SIMD_VEC];
      ALIGNED_float float dat5[2][SIMD_VEC];
      ALIGNED_float float dat6[2][SIMD_VEC];
      ALIGNED_float float dat7[2][SIMD_VEC];
      #pragma simd
      for ( int i=0; i<SIMD_VEC; i++ ){
        int idx1_4 = 4 * arg1.map_data[(n+i) * arg1.map->dim + 0];
        int idx2_4 = 4 * arg1.map_data[(n+i) * arg1.map->dim + 1];
        int idx3_4 = 4 * arg1.map_data[(n+i) * arg1.map->dim + 2];
        int idx5_2 = 2 * arg1.map_data[(n+i) * arg1.map->dim + 0];
        int idx6_2 = 2 * arg1.map_data[(n+i) * arg1.map->dim + 1];
        int idx7_2 = 2 * arg1.map_data[(n+i) * arg1.map->dim + 2];

        dat1[0][i] = (ptr1)[idx1_4 + 0];
        dat1[1][i] = (ptr1)[idx1_4 + 1];
        dat1[2][i] = (ptr1)[idx1_4 + 2];
        dat1[3][i] = (ptr1)[idx1_4 + 3];

        dat2[0][i] = (ptr2)[idx2_4 + 0];
        dat2[1][i] = (ptr2)[idx2_4 + 1];
        dat2[2][i] = (ptr2)[idx2_4 + 2];
        dat2[3][i] = (ptr2)[idx2_4 + 3];

        dat3[0][i] = (ptr3)[idx3_4 + 0];
        dat3[1][i] = (ptr3)[idx3_4 + 1];
        dat3[2][i] = (ptr3)[idx3_4 + 2];
        dat3[3][i] = (ptr3)[idx3_4 + 3];

        dat5[0][i] = (ptr5)[idx5_2 + 0];
        dat5[1][i] = (ptr5)[idx5_2 + 1];

        dat6[0][i] = (ptr6)[idx6_2 + 0];
        dat6[1][i] = (ptr6)[idx6_2 + 1];

        dat7[0][i] = (ptr7)[idx7_2 + 0];
        dat7[1][i] = (ptr7)[idx7_2 + 1];

      }
      #pragma simd
      for ( int i=0; i<SIMD_VEC; i++ ){
        computeGradient_vec(
          &(ptr0)[4 * (n+i)],
          dat1,
          dat2,
          dat3,
          &(ptr4)[2 * (n+i)],
          dat5,
          dat6,
          dat7,
          &(ptr8)[8 * (n+i)],
          &(ptr9)[8 * (n+i)],
          i);
      }
      for ( int i=0; i<SIMD_VEC; i++ ){

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
      int map1idx = arg1.map_data[n * arg1.map->dim + 0];
      int map2idx = arg1.map_data[n * arg1.map->dim + 1];
      int map3idx = arg1.map_data[n * arg1.map->dim + 2];

      computeGradient(
        &(ptr0)[4 * n],
        &(ptr1)[4 * map1idx],
        &(ptr2)[4 * map2idx],
        &(ptr3)[4 * map3idx],
        &(ptr4)[2 * n],
        &(ptr5)[2 * map1idx],
        &(ptr6)[2 * map2idx],
        &(ptr7)[2 * map3idx],
        &(ptr8)[8 * n],
        &(ptr9)[8 * n]);
    }
  }

  if (exec_size == 0 || exec_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[5].name      = name;
  OP_kernels[5].count    += 1;
  OP_kernels[5].time     += wall_t2 - wall_t1;
  OP_kernels[5].transfer += (float)set->size * arg1.size;
  OP_kernels[5].transfer += (float)set->size * arg5.size;
  OP_kernels[5].transfer += (float)set->size * arg0.size;
  OP_kernels[5].transfer += (float)set->size * arg4.size;
  OP_kernels[5].transfer += (float)set->size * arg8.size;
  OP_kernels[5].transfer += (float)set->size * arg9.size;
  OP_kernels[5].transfer += (float)set->size * arg1.map->dim * 4.0f;
}
//
// auto-generated by op2.py
//

//user function
inline void initGaussianLandslide(const float *center, float *values, const float *mesh_xmin, const float *A, const double *t, const float *lx, const float *ly, const float *v) {
  float x = center[0];
  float y = center[1];
  values[3] = (*mesh_xmin-x)*(x<0.0)-5.0*(x>=0.0)+
      *A*(*t<1.0/(*v))*exp(-1.0* *lx* *lx*(x+3.0-*v**t)*(x+3.0-*v**t)-*ly**ly*y*y)
      +*A*(*t>=1.0/(*v))*exp(-*lx*(x+3.0-1.0)**lx*(x+3.0-1.0)-*ly**ly*y*y);
}

// host stub function
void op_par_loop_initGaussianLandslide(char const *name, op_set set,
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
  ALIGNED_float const float * __restrict__ ptr0 = (float *) arg0.data;
  __assume_aligned(ptr0,float_ALIGN);
  ALIGNED_float       float * __restrict__ ptr1 = (float *) arg1.data;
  __assume_aligned(ptr1,float_ALIGN);

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(25);
  op_timers_core(&cpu_t1, &wall_t1);


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  initGaussianLandslide");
  }

  int exec_size = op_mpi_halo_exchanges(set, nargs, args);

  if (exec_size >0) {

    #ifdef VECTORIZE
    #pragma novector
    for ( int n=0; n<(exec_size/SIMD_VEC)*SIMD_VEC; n+=SIMD_VEC ){
      #pragma simd
      for ( int i=0; i<SIMD_VEC; i++ ){
        initGaussianLandslide(
          &(ptr0)[2 * (n+i)],
          &(ptr1)[4 * (n+i)],
          (float*)arg2.data,
          (float*)arg3.data,
          (double*)arg4.data,
          (float*)arg5.data,
          (float*)arg6.data,
          (float*)arg7.data);
      }
      for ( int i=0; i<SIMD_VEC; i++ ){
      }
      for ( int i=0; i<SIMD_VEC; i++ ){
      }
      for ( int i=0; i<SIMD_VEC; i++ ){
      }
      for ( int i=0; i<SIMD_VEC; i++ ){
      }
      for ( int i=0; i<SIMD_VEC; i++ ){
      }
      for ( int i=0; i<SIMD_VEC; i++ ){
      }
    }
    //remainder
    for ( int n=(exec_size/SIMD_VEC)*SIMD_VEC; n<exec_size; n++ ){
    #else
    for ( int n=0; n<exec_size; n++ ){
    #endif
      initGaussianLandslide(
        &(ptr0)[2*n],
        &(ptr1)[4*n],
        (float*)arg2.data,
        (float*)arg3.data,
        (double*)arg4.data,
        (float*)arg5.data,
        (float*)arg6.data,
        (float*)arg7.data);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[25].name      = name;
  OP_kernels[25].count    += 1;
  OP_kernels[25].time     += wall_t2 - wall_t1;
  OP_kernels[25].transfer += (float)set->size * arg0.size;
  OP_kernels[25].transfer += (float)set->size * arg1.size * 2.0f;
}

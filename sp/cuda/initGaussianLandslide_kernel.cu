//
// auto-generated by op2.py
//

//user function
__device__ void initGaussianLandslide_gpu( const float *center, float *values, const float *mesh_xmin, const float *A, const double *t, const float *lx, const float *ly, const float *v) {
  float x = center[0];
  float y = center[1];
  values[3] = (*mesh_xmin-x)*(x<0.0)-5.0*(x>=0.0)+
      *A*(*t<1.0/(*v))*exp(-1.0* *lx* *lx*(x+3.0-*v**t)*(x+3.0-*v**t)-*ly**ly*y*y)
      +*A*(*t>=1.0/(*v))*exp(-*lx*(x+3.0-1.0)**lx*(x+3.0-1.0)-*ly**ly*y*y);

}

// CUDA kernel function
__global__ void op_cuda_initGaussianLandslide(
  const float *__restrict arg0,
  float *arg1,
  const float *arg2,
  const float *arg3,
  const double *arg4,
  const float *arg5,
  const float *arg6,
  const float *arg7,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    initGaussianLandslide_gpu(arg0+n*2,
                          arg1+n*4,
                          arg2,
                          arg3,
                          arg4,
                          arg5,
                          arg6,
                          arg7);
  }
}


//host stub function
void op_par_loop_initGaussianLandslide(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7){

  float*arg2h = (float *)arg2.data;
  float*arg3h = (float *)arg3.data;
  double*arg4h = (double *)arg4.data;
  float*arg5h = (float *)arg5.data;
  float*arg6h = (float *)arg6.data;
  float*arg7h = (float *)arg7.data;
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

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(15);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[15].name      = name;
  OP_kernels[15].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  initGaussianLandslide");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    //transfer constants to GPU
    int consts_bytes = 0;
    consts_bytes += ROUND_UP(1*sizeof(float));
    consts_bytes += ROUND_UP(1*sizeof(float));
    consts_bytes += ROUND_UP(1*sizeof(double));
    consts_bytes += ROUND_UP(1*sizeof(float));
    consts_bytes += ROUND_UP(1*sizeof(float));
    consts_bytes += ROUND_UP(1*sizeof(float));
    reallocConstArrays(consts_bytes);
    consts_bytes = 0;
    arg2.data   = OP_consts_h + consts_bytes;
    arg2.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((float *)arg2.data)[d] = arg2h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(float));
    arg3.data   = OP_consts_h + consts_bytes;
    arg3.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((float *)arg3.data)[d] = arg3h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(float));
    arg4.data   = OP_consts_h + consts_bytes;
    arg4.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((double *)arg4.data)[d] = arg4h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(double));
    arg5.data   = OP_consts_h + consts_bytes;
    arg5.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((float *)arg5.data)[d] = arg5h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(float));
    arg6.data   = OP_consts_h + consts_bytes;
    arg6.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((float *)arg6.data)[d] = arg6h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(float));
    arg7.data   = OP_consts_h + consts_bytes;
    arg7.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((float *)arg7.data)[d] = arg7h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(float));
    mvConstArraysToDevice(consts_bytes);

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_15
      int nthread = OP_BLOCK_SIZE_15;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_initGaussianLandslide<<<nblocks,nthread>>>(
      (float *) arg0.data_d,
      (float *) arg1.data_d,
      (float *) arg2.data_d,
      (float *) arg3.data_d,
      (double *) arg4.data_d,
      (float *) arg5.data_d,
      (float *) arg6.data_d,
      (float *) arg7.data_d,
      set->size );
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[15].time     += wall_t2 - wall_t1;
  OP_kernels[15].transfer += (float)set->size * arg0.size;
  OP_kernels[15].transfer += (float)set->size * arg1.size * 2.0f;
}

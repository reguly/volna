//
// auto-generated by op2.py
//

//user function
__device__ void limiter_gpu( const float *q, float *lim,
                    const float *value, const float *gradient,
                    const float *edgecenter1, const float *edgecenter2,
                    const float *edgecenter3,
                    float *zeroInit,
                    const float *cellcenter) {

  float facevalue[3], dx[3], dy[3];
  int i, j;
  float max[3], edgealpha[3];
  dx[0] = (edgecenter1[0] - cellcenter[0]);
  dy[0] = (edgecenter1[1] - cellcenter[1]);
  dx[1] = (edgecenter2[0] - cellcenter[0]);
  dy[1] = (edgecenter2[1] - cellcenter[1]);
  dx[2] = (edgecenter3[0] - cellcenter[0]);
  dy[2] = (edgecenter3[1] - cellcenter[1]);

  if(q[0] > EPS_cuda){






  for(j=0;j<4;j++){
   for(i =0 ; i<3; i++){
    facevalue[i] = value[j] + (((gradient[2*j]*dx[i]) + (gradient[2*j + 1]*dy[i])));
     if(facevalue[i] > q[2*j + 1]) {
      edgealpha[i] = (q[2*j + 1] - value[j]) / (facevalue[i] - value[j]);
     } else if (facevalue[i] < q[2*j]){
      edgealpha[i] = (q[2*j] - value[j]) / (facevalue[i] - value[j]);
     } else{
      edgealpha[i] = 1.0f;
     }
    max[i] = edgealpha[i] < 1.0f ? edgealpha[i] : 1.0f;
   }
   lim[j] = max[0] < max[1] ? max[0] : max[1];
   lim[j] = lim[j] < max[2] ? lim[j]: max[2];
  }
  lim[0] = lim[0] < lim[1] ? lim[0]: lim[1];
  lim[0] = lim[0] < lim[2] ? lim[0]: lim[2];
  lim[0] = lim[0] < lim[3] ? lim[0]: lim[3];
  } else {
    lim[0] = 0.0f;
    lim[1] = 0.0f;
    lim[2] = 0.0f;
    lim[3] = 0.0f;
  }
  zeroInit[0] = 0.0f;
  zeroInit[1] = 0.0f;
  zeroInit[2] = 0.0f;
  zeroInit[3] = 0.0f;

}

// CUDA kernel function
__global__ void op_cuda_limiter(
  const float *__restrict ind_arg0,
  const int *__restrict opDat4Map,
  const float *__restrict arg0,
  float *arg1,
  const float *__restrict arg2,
  const float *__restrict arg3,
  float *arg7,
  const float *__restrict arg8,
  int start,
  int end,
  int   set_size) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    int map4idx;
    int map5idx;
    int map6idx;
    map4idx = opDat4Map[n + set_size * 0];
    map5idx = opDat4Map[n + set_size * 1];
    map6idx = opDat4Map[n + set_size * 2];

    //user-supplied kernel call
    limiter_gpu(arg0+n*8,
            arg1+n*4,
            arg2+n*4,
            arg3+n*8,
            ind_arg0+map4idx*2,
            ind_arg0+map5idx*2,
            ind_arg0+map6idx*2,
            arg7+n*4,
            arg8+n*2);
  }
}


//host stub function
void op_par_loop_limiter(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8){

  int nargs = 9;
  op_arg args[9];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;
  args[7] = arg7;
  args[8] = arg8;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(23);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[23].name      = name;
  OP_kernels[23].count    += 1;


  int    ninds   = 1;
  int    inds[9] = {-1,-1,-1,-1,0,0,0,-1,-1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: limiter\n");
  }
  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_23
      int nthread = OP_BLOCK_SIZE_23;
    #else
      int nthread = OP_block_size;
    #endif

    for ( int round=0; round<2; round++ ){
      if (round==1) {
        op_mpi_wait_all_grouped(nargs, args, 2);
      }
      int start = round==0 ? 0 : set->core_size;
      int end = round==0 ? set->core_size : set->size + set->exec_size;
      if (end-start>0) {
        int nblocks = (end-start-1)/nthread+1;
        op_cuda_limiter<<<nblocks,nthread>>>(
        (float *)arg4.data_d,
        arg4.map_data_d,
        (float*)arg0.data_d,
        (float*)arg1.data_d,
        (float*)arg2.data_d,
        (float*)arg3.data_d,
        (float*)arg7.data_d,
        (float*)arg8.data_d,
        start,end,set->size+set->exec_size);
      }
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[23].time     += wall_t2 - wall_t1;
}

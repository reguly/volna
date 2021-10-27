//
// auto-generated by op2.py
//

//user function
__device__ void computeFluxes_gpu( const float *cellLeft, const float *cellRight,
                                const float *alphaleft, const float *alpharight,
                                const float *edgeLength, const float *edgeNormals,
                                const float *leftcellCenters, const float *rightcellCenters,
                                const float *edgeCenters,
                                const float *leftGradient, const float *rightGradient,
                                const int *isRightBoundary,
                                float *bathySource, float *out,
                                float *maxEdgeEigenvalues, const float *zmin) {

  float leftCellValues[4];
  float rightCellValues[4];
  float InterfaceBathy;
  float zL, zR;
  float uR,vR,uL,vL;
  leftCellValues[0] = cellLeft[0];
  leftCellValues[1] = cellLeft[1];
  leftCellValues[2] = cellLeft[2];
  leftCellValues[3] = cellLeft[3];
  float dxl, dyl, dxr, dyr;
  dxl = (edgeCenters[0] - leftcellCenters[0]);
  dyl = (edgeCenters[1] - leftcellCenters[1]);
  dxr = (edgeCenters[0] - rightcellCenters[0]);
  dyr = (edgeCenters[1] - rightcellCenters[1]);

  leftCellValues[0] += alphaleft[0] * ((dxl * leftGradient[0])+(dyl * leftGradient[1]));
  leftCellValues[1] += alphaleft[0] * ((dxl * leftGradient[2])+(dyl * leftGradient[3]));
  leftCellValues[2] += alphaleft[0] * ((dxl * leftGradient[4])+(dyl * leftGradient[5]));
  leftCellValues[3] += alphaleft[0] * ((dxl * leftGradient[6])+(dyl * leftGradient[7]));
  if (leftCellValues[0] >= 1e-3){
     uL = leftCellValues[1]/leftCellValues[0];
     vL = leftCellValues[2]/leftCellValues[0];
  } else {
     uL = 0.0f;
     vL = 0.0f;
  }
  zL = cellLeft[3] - cellLeft[0];
  zR = cellRight[3] - cellRight[0];
  if (!*isRightBoundary) {
    rightCellValues[0] = cellRight[0];
    rightCellValues[1] = cellRight[1];
    rightCellValues[2] = cellRight[2];
    rightCellValues[3] = cellRight[3];

    rightCellValues[0] += alpharight[0] * ((dxr * rightGradient[0])+(dyr * rightGradient[1]));
    rightCellValues[1] += alpharight[0] * ((dxr * rightGradient[2])+(dyr * rightGradient[3]));
    rightCellValues[2] += alpharight[0] * ((dxr * rightGradient[4])+(dyr * rightGradient[5]));
    rightCellValues[3] += alpharight[0] * ((dxr * rightGradient[6])+(dyr * rightGradient[7]));
    if (rightCellValues[0] >= 1e-3){
       uR = rightCellValues[1]/rightCellValues[0];
       vR = rightCellValues[2]/rightCellValues[0];
     } else {
       uR = 0.0f;
       vR = 0.0f;
     }
   } else {
    float nx = edgeNormals[0];
    float ny = edgeNormals[1];
    float inNormalVelocity = uL * nx + vL * ny;
    float inTangentVelocity = -1.0f * uL * ny + vL * nx;
    float outNormalVelocity = 0.0f;
    float outTangentVelocity = 0.0f;

    float wet = (fabs(*zmin) - zL) > 0.0f ? (fabs(*zmin) - zL) : EPS_cuda;
    float critical = sqrt(uL*uL + vL*vL);
    outTangentVelocity = inTangentVelocity;
    if (critical < sqrt(g_cuda*leftCellValues[0])){
      rightCellValues[0] = wet;
      rightCellValues[3] = wet + zL;
      outNormalVelocity = 0.0f;
    } else {
      rightCellValues[0] = leftCellValues[0];
      rightCellValues[3] = leftCellValues[3];
      outNormalVelocity = inNormalVelocity;
    }

    uR = outNormalVelocity * nx - outTangentVelocity * ny;
    vR = outNormalVelocity * ny + outTangentVelocity * nx;
    rightCellValues[1] = uR*rightCellValues[0];
    rightCellValues[2] = vR*rightCellValues[0];
    zR = zL;
  }

  rightCellValues[3] -= rightCellValues[0];
  leftCellValues[3] -= leftCellValues[0];



  InterfaceBathy = leftCellValues[3] > rightCellValues[3] ? leftCellValues[3] : rightCellValues[3];
  bathySource[0] =0.5f * g_cuda * (leftCellValues[0]*leftCellValues[0]);
  bathySource[1] =0.5f * g_cuda * (rightCellValues[0]*rightCellValues[0]);
  float hL = (leftCellValues[0] + leftCellValues[3] - InterfaceBathy);
  hL = hL > 0.0f ? hL : 0.0f;
  float hR = (rightCellValues[0] + rightCellValues[3] - InterfaceBathy);
  hR = hR > 0.0f ? hR : 0.0f;
  bathySource[0] -= .5f * g_cuda * (hL * hL);
  bathySource[1] -= .5f * g_cuda * (hR * hR);

  bathySource[2] = -.5f * g_cuda *(leftCellValues[0] + cellLeft[0])*(leftCellValues[3] - zL);
  bathySource[3] = -.5f * g_cuda *(rightCellValues[0] + cellRight[0])*(rightCellValues[3] - zR);

  bathySource[0] *= *edgeLength;
  bathySource[1] *= *edgeLength;
  bathySource[2] *= *edgeLength;
  bathySource[3] *= *edgeLength;



  float cL = sqrt(g_cuda * hL);
  cL = cL > 0.0f ? cL : 0.0f;
  float cR = sqrt(g_cuda * hR);
  cR = cR > 0.0f ? cR : 0.0f;

  float uLn = uL * edgeNormals[0] + vL * edgeNormals[1];
  float uRn = uR * edgeNormals[0] + vR * edgeNormals[1];

  float unStar = 0.5f * (uLn + uRn) + (cL-cR);
  float cStar = 0.5f * (cL + cR) - 0.25f* (uRn-uLn);
  float sL = (uLn - cL) < (unStar - cStar) ? (uLn - cL) : (unStar - cStar);
  float sR = (uRn + cR) > (unStar + cStar) ? (uRn + cR) : (unStar + cStar);

  float sStar;
  sStar = (sL*hR*(uRn - sR) - sR*hL*(uLn - sL))/
          (hR*(uRn - sR) - hL*(uLn - sL));

  if ((leftCellValues[0] <= EPS_cuda) && (rightCellValues[0] > EPS_cuda)) {
      sL = uRn - 2.0f*cR;
      sR = uRn + cR;
      sStar = sL;
  }
  if ((rightCellValues[0] <= EPS_cuda) && (leftCellValues[0] > EPS_cuda)) {
      sR = uLn + 2.0f*cL;
      sL =  uLn - cL;
      sStar = sR;
  }


  float uLp = vL*edgeNormals[0] - uL*edgeNormals[1];
  float uRp = vR*edgeNormals[0] - uR*edgeNormals[1];

  float LeftFluxes_H, LeftFluxes_N, LeftFluxes_U, LeftFluxes_V;
  float HuDotN = (hL*uL) * edgeNormals[0] + (hL*vL) * edgeNormals[1];

  LeftFluxes_H = HuDotN;
  LeftFluxes_U = HuDotN * uL;
  LeftFluxes_V = HuDotN * vL;

  LeftFluxes_N = HuDotN * uLn;

  LeftFluxes_U += (.5f * g_cuda * edgeNormals[0] ) * ( hL * hL );
  LeftFluxes_V += (.5f * g_cuda * edgeNormals[1] ) * ( hL * hL );
  LeftFluxes_N += (.5f * g_cuda ) * ( hL * hL );

  float RightFluxes_H,RightFluxes_N, RightFluxes_U, RightFluxes_V;
  HuDotN = (hR*uR) * edgeNormals[0] + (hR*vR) * edgeNormals[1];

  RightFluxes_H =   HuDotN;
  RightFluxes_U =   HuDotN * uR;
  RightFluxes_V =   HuDotN * vR;

  RightFluxes_N =   HuDotN * uRn;

  RightFluxes_U += (.5f * g_cuda * edgeNormals[0] ) * ( hR * hR );
  RightFluxes_V += (.5f * g_cuda * edgeNormals[1] ) * ( hR * hR );
  RightFluxes_N += (.5f * g_cuda ) * ( hR * hR );














  float sLMinus = sL < 0.0f ? sL : 0.0f;
  float sRPlus = sR > 0.0f ? sR : 0.0f;
  float sRMinussL = sRPlus - sLMinus;
  sRMinussL = sRMinussL < EPS_cuda ?  EPS_cuda : sRMinussL;
  float t1 = sRPlus / sRMinussL;
  float t2 = ( -1.0 * sLMinus ) / sRMinussL;
  float t3 = ( sRPlus * sLMinus ) / sRMinussL;
  float FStar[3];
  FStar[0] =
  ( t1 * LeftFluxes_H ) +
  ( t2 * RightFluxes_H ) +
  ( t3 * ( hR - hL ) );

  FStar[1] =
  ( t1 * LeftFluxes_N ) +
  ( t2 * RightFluxes_N ) +
  ( t3 * ( (hR * uRn) -
          (hL * uLn) ) );

  if( sL >= 0.0f) {
    out[0] = t1*LeftFluxes_H;
    out[1] = t1*LeftFluxes_U;
    out[2] = t1*LeftFluxes_V;
  } else if ((sL < 0.0f) && (sStar >= 0.0f)){
    out[0] = FStar[0];
    FStar[2] = FStar[0] * uLp;
    out[1] = FStar[1]*edgeNormals[0] - FStar[2]*edgeNormals[1];
    out[2] = FStar[1]*edgeNormals[1] + FStar[2]*edgeNormals[0];
  } else if((sStar < 0.0f) && (sR >= 0.0f)){
    out[0] = FStar[0];
    FStar[2] = FStar[0] * uRp;
    out[1] = FStar[1]*edgeNormals[0] - FStar[2]*edgeNormals[1];
    out[2] = FStar[1]*edgeNormals[1] + FStar[2]*edgeNormals[0];
  } else {
    out[0] = t2*RightFluxes_H;
    out[1] = t2*RightFluxes_U;
    out[2] = t2*RightFluxes_V;
  }
  out[0] *= *edgeLength;
  out[1] *= *edgeLength;
  out[2] *= *edgeLength;

  float maximum = fabs(uLn + cL);
  maximum = maximum > fabs(uLn - cL) ? maximum : fabs(uLn - cL);
  maximum = maximum > fabs(uRn + cR) ? maximum : fabs(uRn + cR);
  maximum = maximum > fabs(uRn - cR) ? maximum : fabs(uRn - cR);
  *maxEdgeEigenvalues = maximum;

}

// CUDA kernel function
__global__ void op_cuda_computeFluxes(
  const float *__restrict ind_arg0,
  const float *__restrict ind_arg1,
  const float *__restrict ind_arg2,
  const float *__restrict ind_arg3,
  const int *__restrict opDat0Map,
  const float *__restrict arg4,
  const float *__restrict arg5,
  const float *__restrict arg8,
  const int *__restrict arg11,
  float *arg12,
  float *arg13,
  float *arg14,
  const float *arg15,
  int    block_offset,
  int   *blkmap,
  int   *offset,
  int   *nelems,
  int   *ncolors,
  int   *colors,
  int   nblocks,
  int   set_size) {


  __shared__ int    nelem, offset_b;

  extern __shared__ char shared[];

  if (blockIdx.x+blockIdx.y*gridDim.x >= nblocks) {
    return;
  }
  if (threadIdx.x==0) {

    //get sizes and shift pointers and direct-mapped data

    int blockId = blkmap[blockIdx.x + blockIdx.y*gridDim.x  + block_offset];

    nelem    = nelems[blockId];
    offset_b = offset[blockId];

  }
  __syncthreads(); // make sure all of above completed

  for ( int n=threadIdx.x; n<nelem; n+=blockDim.x ){
    int map0idx;
    int map1idx;
    map0idx = opDat0Map[n + offset_b + set_size * 0];
    map1idx = opDat0Map[n + offset_b + set_size * 1];


    //user-supplied kernel call
    computeFluxes_gpu(ind_arg0+map0idx*4,
                  ind_arg0+map1idx*4,
                  ind_arg1+map0idx*4,
                  ind_arg1+map1idx*4,
                  arg4+(n+offset_b)*1,
                  arg5+(n+offset_b)*2,
                  ind_arg2+map0idx*2,
                  ind_arg2+map1idx*2,
                  arg8+(n+offset_b)*2,
                  ind_arg3+map0idx*8,
                  ind_arg3+map1idx*8,
                  arg11+(n+offset_b)*1,
                  arg12+(n+offset_b)*4,
                  arg13+(n+offset_b)*3,
                  arg14+(n+offset_b)*1,
                  arg15);
  }
}


//host stub function
void op_par_loop_computeFluxes(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9,
  op_arg arg10,
  op_arg arg11,
  op_arg arg12,
  op_arg arg13,
  op_arg arg14,
  op_arg arg15){

  float*arg15h = (float *)arg15.data;
  int nargs = 16;
  op_arg args[16];

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
  args[10] = arg10;
  args[11] = arg11;
  args[12] = arg12;
  args[13] = arg13;
  args[14] = arg14;
  args[15] = arg15;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(6);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[6].name      = name;
  OP_kernels[6].count    += 1;


  int    ninds   = 4;
  int    inds[16] = {0,0,1,1,-1,-1,2,2,-1,3,3,-1,-1,-1,-1,-1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: computeFluxes\n");
  }

  //get plan
  #ifdef OP_PART_SIZE_6
    int part_size = OP_PART_SIZE_6;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    op_plan *Plan = op_plan_get(name,set,part_size,nargs,args,ninds,inds);

    //transfer constants to GPU
    int consts_bytes = 0;
    consts_bytes += ROUND_UP(1*sizeof(float));
    reallocConstArrays(consts_bytes);
    consts_bytes = 0;
    arg15.data   = OP_consts_h + consts_bytes;
    arg15.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((float *)arg15.data)[d] = arg15h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(float));
    mvConstArraysToDevice(consts_bytes);

    //execute plan

    int block_offset = 0;
    for ( int col=0; col<Plan->ncolors; col++ ){
      if (col==Plan->ncolors_core) {
        op_mpi_wait_all_grouped(nargs, args, 2);
      }
      #ifdef OP_BLOCK_SIZE_6
      int nthread = OP_BLOCK_SIZE_6;
      #else
      int nthread = OP_block_size;
      #endif

      dim3 nblocks = dim3(Plan->ncolblk[col] >= (1<<16) ? 65535 : Plan->ncolblk[col],
      Plan->ncolblk[col] >= (1<<16) ? (Plan->ncolblk[col]-1)/65535+1: 1, 1);
      if (Plan->ncolblk[col] > 0) {
        op_cuda_computeFluxes<<<nblocks,nthread>>>(
        (float *)arg0.data_d,
        (float *)arg2.data_d,
        (float *)arg6.data_d,
        (float *)arg9.data_d,
        arg0.map_data_d,
        (float*)arg4.data_d,
        (float*)arg5.data_d,
        (float*)arg8.data_d,
        (int*)arg11.data_d,
        (float*)arg12.data_d,
        (float*)arg13.data_d,
        (float*)arg14.data_d,
        (float*)arg15.data_d,
        block_offset,
        Plan->blkmap,
        Plan->offset,
        Plan->nelems,
        Plan->nthrcol,
        Plan->thrcol,
        Plan->ncolblk[col],
        set->size+set->exec_size);

      }
      block_offset += Plan->ncolblk[col];
    }
    OP_kernels[6].transfer  += Plan->transfer;
    OP_kernels[6].transfer2 += Plan->transfer2;
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  if (OP_diags>1) {
    cutilSafeCall(cudaDeviceSynchronize());
  }
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[6].time     += wall_t2 - wall_t1;
}

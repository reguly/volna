//
// auto-generated by op2.py
//

//user function
inline void gatherLocations(const float *values, float *dest) {
	dest[0] = values[0] + values[3];
  dest[1] = values[0];
  dest[2] = values[1];
  dest[3] = values[2];
  dest[4] = values[3];
}
#ifdef VECTORIZE
//user function -- modified for vectorisation
void gatherLocations_vec( const float values[*][SIMD_VEC], float *dest, int idx ) {
	dest[0] = values[0][idx] + values[3][idx];
  dest[1] = values[0][idx];
  dest[2] = values[1][idx];
  dest[3] = values[2][idx];
  dest[4] = values[3][idx];
}
#endif

// host stub function
void op_par_loop_gatherLocations(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1){

  int nargs = 2;
  op_arg args[2];

  args[0] = arg0;
  args[1] = arg1;
  //create aligned pointers for dats
  ALIGNED_float const float * __restrict__ ptr0 = (float *) arg0.data;
  __assume_aligned(ptr0,float_ALIGN);
  ALIGNED_float       float * __restrict__ ptr1 = (float *) arg1.data;
  __assume_aligned(ptr1,float_ALIGN);

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(9);
  op_timers_core(&cpu_t1, &wall_t1);

  int  ninds   = 1;
  int  inds[2] = {0,-1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: gatherLocations\n");
  }

  #ifdef OP_PART_SIZE_9
    int part_size = OP_PART_SIZE_9;
  #else
    int part_size = OP_part_size;
  #endif


  int set_size = op_mpi_halo_exchanges(set, nargs, args);

  if (set_size >0) {

    // get plan
    op_plan *Plan = op_plan_get_stage_upload(name,set,part_size,nargs,args,ninds,inds,OP_STAGE_ALL,0);

    // execute plan
    int block_offset = 0;
    for ( int col=0; col<Plan->ncolors; col++ ){
      if (col==Plan->ncolors_core) {
        op_mpi_wait_all(nargs, args);
      }
      int nblocks = Plan->ncolblk[col];

      #pragma omp parallel for
      for ( int blockIdx=0; blockIdx<nblocks; blockIdx++ ){
        int blockId  = Plan->blkmap[blockIdx + block_offset];
        int nelem    = Plan->nelems[blockId];
        int offset_b = Plan->offset[blockId];
        #ifdef VECTORIZE
        //peel left remainder
        for ( int n=offset_b; n<((offset_b-1)/SIMD_VEC+1)*SIMD_VEC; n++ ){
          int map0idx = arg0.map_data[n * arg0.map->dim + 0];

          gatherLocations(
            &(ptr0)[4 * map0idx],
            &(ptr1)[5 * n]);
        }
        #pragma novector
        for ( int n=((offset_b-1)/SIMD_VEC+1)*SIMD_VEC; n<((offset_b+nelem)/SIMD_VEC)*SIMD_VEC; n+=SIMD_VEC ){
          if (n+SIMD_VEC >= set->core_size) {
            op_mpi_wait_all(nargs, args);
          }
          ALIGNED_float float dat0[4][SIMD_VEC];
          #pragma omp simd aligned(ptr0,ptr1)
          for ( int i=0; i<SIMD_VEC; i++ ){
            int idx0_4 = 4 * arg0.map_data[(n+i) * arg0.map->dim + 0];

            dat0[0][i] = (ptr0)[idx0_4 + 0];
            dat0[1][i] = (ptr0)[idx0_4 + 1];
            dat0[2][i] = (ptr0)[idx0_4 + 2];
            dat0[3][i] = (ptr0)[idx0_4 + 3];

          }
          #pragma omp simd aligned(ptr0,ptr1)
          for ( int i=0; i<SIMD_VEC; i++ ){
            gatherLocations_vec(
              dat0,
              &(ptr1)[5 * (n+i)],
              i);
          }
          for ( int i=0; i<SIMD_VEC; i++ ){

          }
        }

        //remainder
        for ( int n=((offset_b+nelem)/SIMD_VEC)*SIMD_VEC; n<offset_b+nelem; n++ ){
        #else
        #pragma omp simd aligned(ptr0,ptr1)
        for ( int n=offset_b; n<offset_b+nelem; n++ ){
        #endif
          int map0idx = arg0.map_data[n * arg0.map->dim + 0];

          gatherLocations(
            &(ptr0)[4 * map0idx],
            &(ptr1)[5 * n]);
        }
      }
      block_offset += nblocks;
    }
  }

  if (set_size == 0 || set_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[9].name      = name;
  OP_kernels[9].count    += 1;
  OP_kernels[9].time     += wall_t2 - wall_t1;
  OP_kernels[9].transfer += (float)set->size * arg0.size;
  OP_kernels[9].transfer += (float)set->size * arg1.size;
  OP_kernels[9].transfer += (float)set->size * arg0.map->dim * 4.0f;
}
#undef VECTORIZE

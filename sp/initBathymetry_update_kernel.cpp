//
// auto-generated by op2.py on 2014-10-07 13:54
//

//user function
#include "initBathymetry_update.h"

// host stub function
void op_par_loop_initBathymetry_update(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1){

  int nargs = 2;
  op_arg args[2];

  args[0] = arg0;
  args[1] = arg1;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(11);
  op_timers_core(&cpu_t1, &wall_t1);


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  initBathymetry_update");
  }

  op_mpi_halo_exchanges(set, nargs, args);
  // set number of threads
  #ifdef _OPENMP
    int nthreads = omp_get_max_threads();
  #else
    int nthreads = 1;
  #endif

  if (set->size >0) {

    // execute plan
    #pragma omp parallel for
    for ( int thr=0; thr<nthreads; thr++ ){
      int start  = (set->size* thr)/nthreads;
      int finish = (set->size*(thr+1))/nthreads;
      for ( int n=start; n<finish; n++ ){
        initBathymetry_update(
          &((float*)arg0.data)[4*n],
          (int*)arg1.data);
      }
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[11].name      = name;
  OP_kernels[11].count    += 1;
  OP_kernels[11].time     += wall_t2 - wall_t1;
  OP_kernels[11].transfer += (float)set->size * arg0.size * 2.0f;
}

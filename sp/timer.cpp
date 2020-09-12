#include "timer.h"
#include <sys/time.h>

int sl_kern_max = 0;
sl_kernel* sl_kernels;
int sl_kernel_max = 0;

int sl_is_root() { return 1; }

void* sl_realloc(void *ptr, size_t size) {
  return realloc(ptr, size);
}

void sl_timing_realloc_manytime(int kernel, int num_timers) {

  
  int sl_kern_max_new;

  if (kernel >= sl_kern_max) {
    sl_kern_max_new = kernel + 10;
    sl_kernels = (sl_kernel *)sl_realloc(sl_kernels,
                                         sl_kern_max_new * sizeof(sl_kernel));
    if (sl_kernels == NULL) {
      printf(" op_timing_realloc error \n");
      exit(-1);
    }
    printf("kernel=%d, num=%d\n", kernel, num_timers);
    for (int n = sl_kern_max; n < sl_kern_max_new; n++) {

      sl_kernels[n].counts = (int*)malloc(num_timers * sizeof(int));
      sl_kernels[n].times = (double*)malloc(num_timers * sizeof(double));
      sl_kernels[n].cpu_t1 = (double*)malloc(num_timers * sizeof(double));
      sl_kernels[n].cpu_t2 = (double*)malloc(num_timers * sizeof(double));
      sl_kernels[n].wall_t1 = (double*)malloc(num_timers * sizeof(double));
      sl_kernels[n].wall_t2 = (double*)malloc(num_timers * sizeof(double));
      for (int t = 0; t < num_timers; t++) {

        sl_kernels[n].counts[t] = 0;
        sl_kernels[n].times[t] = 0.0f;
        sl_kernels[n].cpu_t1[t] = 0.0f;
        sl_kernels[n].cpu_t2[t] = 0.0f;
        sl_kernels[n].wall_t1[t] = 0.0f;
        sl_kernels[n].wall_t2[t] = 0.0f;
      }
      sl_kernels[n].ntimes = num_timers;
      sl_kernels[n].name = "unused";
    }
    sl_kern_max = sl_kern_max_new;
  }
}

void sl_timers_core(double *cpu, double *et) {
  (void)cpu;
  struct timeval t;

  gettimeofday(&t, (struct timezone *)0);
  *et = t.tv_sec + t.tv_usec * 1.0e-6;
}



void sl_timing_output_core() {
  printf("\n "
          "-----------------------------------------------------------------"
          "--------------------------\n");
  printf("\n  thread id  count    threadtime(std)    kernel name\n\n");
  for (int n = 0; n < sl_kern_max; n++) {
    if(strcmp(sl_kernels[n].name, "unused") != 0){
      if (sl_is_root()){
        for (int i=0; i<sl_kernels[n].ntimes; i++) {

          if(i == 0){
            printf(" %2d       %6d;      %8.4f;       %s\n", i, sl_kernels[n].counts[i], sl_kernels[n].times[i], 
            sl_kernels[n].name);
          }else{
            printf(" %2d       %6d;      %8.4f;\n", i, sl_kernels[n].counts[i], sl_kernels[n].times[i]);
          }

        }
        printf("\n\n");
      }
    }

  }
}
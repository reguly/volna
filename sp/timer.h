#ifndef __SL_TIMER_H
#define __SL_TIMER_H

#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/queue.h> //contains double linked list implementation

typedef struct {
  char const *name;
  int ntimes;
  int* counts;
  double* times;
  double* cpu_t1;
  double* cpu_t2;
  double* wall_t1;
  double* wall_t2;
} sl_kernel;


void* sl_realloc(void *ptr, size_t size);
void sl_timing_realloc_manytime(int kernel, int num_timers);
void sl_timers_core(double *cpu, double *et);
void sl_timing_output_core();

#endif
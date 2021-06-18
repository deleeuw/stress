#include <stdio.h>
#include <omp.h>

void hellomp() {
  int tid, nthreads;
#pragma omp parallel private(tid)
  {
    nthreads = omp_get_num_threads();
    tid = omp_get_thread_num();
    printf("Hello from thread %d of %d\n", tid, nthreads);
  }
}
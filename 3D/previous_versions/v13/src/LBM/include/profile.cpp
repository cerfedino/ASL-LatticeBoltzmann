#include "profile.h"
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef PROFILE
profiler *init_profiler(uint64_t flops, uint64_t bytes) { return NULL; }
void start_run(profiler *p) {}
void end_run(profiler *p) {}
profiler_stats finish_profiler(profiler *p) {
  profiler_stats tmp;
  return tmp;
}
#else
// Returns the current cycle time as accurately as possible
inline uint64_t time() {
  uint64_t t;
  __asm__ volatile("cpuid;"
                   "mfence;"
                   "rdtsc;"
                   "shl  $32, %%rdx;"
                   "lea  (%%rax, %%rdx),  %0;"
                   : "=r"(t)
                   :
                   : "rax", "rdx", "rbx", "rcx", "memory");
  return t;
}

profiler *init_profiler(uint64_t flops, uint64_t bytes) {
  profiler *p = (profiler *)malloc(sizeof(profiler));
  p->_cycles = 0;
  p->_runs = 0;
  p->_min_cycles = UINT64_MAX;
  p->flops = flops;
  p->bytes = bytes;
  return p;
}

void start_run(profiler *p) {
  p->_runs++;
  p->_start = time();
}

void end_run(profiler *p) {
  uint64_t end_t = time();
  uint64_t time_diff = end_t - p->_start;
  if (time_diff < p->_min_cycles) {
    p->_min_cycles = time_diff;
  }
  if (end_t < p->_start) {
    printf("End time is before than start time, disregarding run... end_t: %llu, p->_start: %llu, time: %llu\n", end_t, p->_start, time());
    p->_runs--;
    return;
  }
  if (p->_cycles > (UINT64_MAX - time_diff)) {
    printf("\nOVERFLOW WHEN ENDING PROFILE RUN, time_diff: %llu, end_t: %llu, p->_start: %llu, time: %llu\n", time_diff, end_t, p->_start, time());
    exit(1);
  }
  p->_cycles += end_t - p->_start;
}

profiler_stats finish_profiler(profiler *p) {
  profiler_stats ps;
  ps.total_cycles = p->_cycles;
  ps.cycles = p->_min_cycles;
  ps.runs = p->_runs;
  ps.performance = (double)(p->flops) / p->_min_cycles;
  ps.arithmetic_intensity = (double)(p->flops) / p->bytes;
  return ps;
}
#endif
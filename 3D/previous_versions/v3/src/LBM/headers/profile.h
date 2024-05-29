#pragma once

#include <stdint.h>
#define PROFILE
typedef struct profiler {
  // Flops that will be profiled
  uint64_t flops;
  // Bytes that will be profiled
  uint64_t bytes;

  // How many runs were performed so far
  int _runs;
  // Start time
  uint64_t _start;
  // How many cycles were counted so far
  uint64_t _cycles;
  // The minimum amount of cycles encountered
  uint64_t _min_cycles;
} profiler;

typedef struct profiler_stats {
  int runs;
  uint64_t cycles;
  uint64_t total_cycles;
  double performance;
  double arithmetic_intensity;
} profiler_stats;

profiler *init_profiler(uint64_t flops, uint64_t bytes);
void start_run(profiler *p);
void end_run(profiler *p);
profiler_stats finish_profiler(profiler *p);
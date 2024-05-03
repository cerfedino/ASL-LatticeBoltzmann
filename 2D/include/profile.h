#pragma once

#include <stdint.h>

typedef struct profiler {
    // Flops that will be profiled
    uint64_t flops;

    // How many runs were performed so far
    int _runs;
    // Start time
    uint64_t _start;
    // How many cycles were counted so far
    uint64_t _cycles;
} profiler;

typedef struct profiler_stats {
    int runs;
    uint64_t cycles;
    double performance;
} profiler_stats;

profiler *init_profiler(uint64_t flops);
void start_run(profiler *p);
void end_run(profiler *p);
profiler_stats finish_profiler(profiler *p);
#include "profile.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef PROFILE
profiler *init_profiler(uint64_t flops) { return NULL; }
void start_run(profiler *p) {}
void end_run(profiler *p) {}
profiler_stats finish_profiler(profiler *p) {
    profiler_stats tmp;
    return tmp;
}
#else
// Returns the current cycle time as accurately as possible
uint64_t time() {
    uint64_t t;
    __asm__ volatile("cpuid;"
                     "lfence;"
                     "rdtsc;"
                     "shl  $32, %%rdx;"
                     "lea  (%%rax, %%rdx),  %0;"
                     : "=r"(t)
                     :
                     : "rax", "rdx");
    return t;
}

profiler *init_profiler(uint64_t flops) {
    profiler *p = (profiler *)malloc(sizeof(profiler));
    p->_cycles = 0;
    p->_runs = 0;
    p->flops = flops;
    return p;
}

void start_run(profiler *p) {
    p->_runs++;
    p->_start = time();
}

void end_run(profiler *p) {
    uint64_t end_t = time();
    p->_cycles += end_t - p->_start;
}

profiler_stats finish_profiler(profiler *p) {
    profiler_stats ps;
    ps.cycles = p->_cycles;
    ps.runs = p->_runs;
    ps.performance = (double)p->_cycles / (p->_runs * p->flops);
    return ps;
}
#endif
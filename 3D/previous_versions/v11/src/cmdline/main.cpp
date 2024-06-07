#include <iostream>
#include <stdio.h>
// Now Linux only.
#include "profile.h"
#include "utils.h"
#include <fstream>
#include <iostream>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#ifdef BENCHMARK
#include <papi.h>
#endif

// the variables are in utils.h
using namespace std;

profiler *compute_density_momentum_profiler, *collision_profiler, *stream_profiler;

int time_lbm = 0;
int time_lbm_x = 0;
int time_lbm_y = 0;
int time_lbm_z = 0;

#ifdef MNx
#define NX (MNx)
#endif
#ifdef MNy
#define NY (MNy)
#endif
#ifdef MNz
#define NZ (MNz)
#endif
#ifdef MNt
#define NT (MNt)
#endif

double *density_field;
double *velocity_field_x;
double *velocity_field_y;
double *velocity_field_z;
double *particle_distributions;
int direction_size = 15;
#define c_s_4 (2 * c_s * c_s * c_s * c_s)
#define c_s_2 (c_s * c_s)
#define c_s_2_2 (2 * c_s * c_s)
#define weights_172 (1.0 / 72.0)
#define weights_19 (1.0 / 9.0)
#define weights_19 (1.0 / 9.0)
#define weights_29 (2.0 / 9.0)

#ifdef BENCHMARK
int papi_event_set = PAPI_NULL;

long long papi_values[2];

long long papi_density_values[2] = {0, 0};
long long papi_collision_values[2] = {0, 0};
long long papi_stream_values[2] = {0, 0};
#endif

inline int scalar_index(int x, int y, int z) { return (z * NX * NY) + (y * NX) + x; }
inline int scalar_index_1(int x, int y, int z, int w) { return ((NX + x - time_lbm_x) % NX + y * NX + z * NX * NY + w * NX * NY * NZ); }
inline int scalar_index_2(int x, int y, int z, int w) { return ((x + time_lbm_x) % NX + y * NX + z * NX * NY + w * NX * NY * NZ); }
inline int scalar_index_3(int x, int y, int z, int w) { return (x + ((NY + y - time_lbm_y) % NY) * NX + z * NX * NY + w * NX * NY * NZ); }
inline int scalar_index_4(int x, int y, int z, int w) { return (x + ((NY + y + time_lbm_y) % NY) * NX + z * NX * NY + w * NX * NY * NZ); }
inline int scalar_index_5(int x, int y, int z, int w) { return (x + y * NX + ((NZ + z - time_lbm_z) % NZ) * NX * NY + w * NX * NY * NZ); }
inline int scalar_index_6(int x, int y, int z, int w) { return (x + y * NX + ((NZ + z + time_lbm_z) % NZ) * NX * NY + w * NX * NY * NZ); }
inline int scalar_index_7(int x, int y, int z, int w) { return (((NX + x - time_lbm_x) % NX) + ((NY + y - time_lbm_y) % NY) * NX + (NZ + z - time_lbm_z) % NZ * NX * NY + w * NX * NY * NZ); }
inline int scalar_index_8(int x, int y, int z, int w) { return (((NX + x + time_lbm_x) % NX) + ((NY + y + time_lbm_y) % NY) * NX + (NZ + z + time_lbm_z) % NZ * NX * NY + w * NX * NY * NZ); }
inline int scalar_index_9(int x, int y, int z, int w) { return (((NX + x - time_lbm_x) % NX) + ((NY + y - time_lbm_y) % NY) * NX + (NZ + z + time_lbm_z) % NZ * NX * NY + w * NX * NY * NZ); }
inline int scalar_index_10(int x, int y, int z, int w) { return (((NX + x + time_lbm_x) % NX) + ((NY + y + time_lbm_y) % NY) * NX + (NZ + z - time_lbm_z) % NZ * NX * NY + w * NX * NY * NZ); }
inline int scalar_index_11(int x, int y, int z, int w) { return (((NX + x - time_lbm_x) % NX) + ((NY + y + time_lbm_y) % NY) * NX + (NZ + z - time_lbm_z) % NZ * NX * NY + w * NX * NY * NZ); }
inline int scalar_index_12(int x, int y, int z, int w) { return (((NX + x + time_lbm_x) % NX) + ((NY + y - time_lbm_y) % NY) * NX + (NZ + z + time_lbm_z) % NZ * NX * NY + w * NX * NY * NZ); }
inline int scalar_index_13(int x, int y, int z, int w) { return (((NX + x + time_lbm_x) % NX) + ((NY + y - time_lbm_y) % NY) * NX + (NZ + z - time_lbm_z) % NZ * NX * NY + w * NX * NY * NZ); }
inline int scalar_index_14(int x, int y, int z, int w) { return (((NX + x - time_lbm_x) % NX) + ((NY + y + time_lbm_y) % NY) * NX + (NZ + z + time_lbm_z) % NZ * NX * NY + w * NX * NY * NZ); }
inline int scalar_index(int x, int y, int z, int w) { return (x + y * NX + z * NX * NY + w * NX * NY * NZ); }

void initialise() {
  int box_flatten_length = NX * NY * NZ;
  int distributions_flatten_length = box_flatten_length * direction_size;
  density_field = (double *)malloc(box_flatten_length * sizeof(double));
  velocity_field_x = (double *)malloc(box_flatten_length * sizeof(double));
  velocity_field_y = (double *)malloc(box_flatten_length * sizeof(double));
  velocity_field_z = (double *)malloc(box_flatten_length * sizeof(double));
  particle_distributions = (double *)malloc(distributions_flatten_length * sizeof(double));
  srand(42);

  for (int i = 0; i < NX * NY * NZ; i++) {
    density_field[i] = 1.0;
  }
  for (int x = 0; x < NX; x++) {
    for (int y = 0; y < NY; y++) {
      for (int z = 0; z < NZ; z++) {
        for (int i = 0; i < direction_size; i++) {
          if (i == 0)
            particle_distributions[scalar_index(x, y, z, i)] = weights_29;
          else if (i < 7)
            particle_distributions[scalar_index(x, y, z, i)] = weights_19;
          else
            particle_distributions[scalar_index(x, y, z, i)] = weights_172;
        }

        velocity_field_x[scalar_index(x, y, z)] = 0.05 * (z % 15);
        velocity_field_y[scalar_index(x, y, z)] = 0.025 * (z % 15);
        velocity_field_z[scalar_index(x, y, z)] = 0.025 * (z % 15);
      }
    }
  }

  compute_density_momentum_profiler = init_profiler(NX * NY * NZ * direction_size * 7 + NX * NY * NZ * 3, 8 * NX * NY * NZ * direction_size * 7 + (NX * NY * NZ * 3 + 15) * 8);
  collision_profiler = init_profiler((NX * NY * NZ * direction_size * 29 + 2), 8 * (4 * NX * NY * NZ + 15 + 15));
  stream_profiler = init_profiler(NX * NZ * 5 * 6, NX * NY * NZ * direction_size * 8 * 2);

#ifdef BENCHMARK
  papi_init(&papi_event_set);
#endif
}

void compute_density_momentum_moment() {

  for (int z = 0; z < NZ; z++) {
    for (int y = 0; y < NY; y++) {
      for (int x = 0; x < NX; x++) {
        double new_density = particle_distributions[scalar_index(x, y, z, 0)] + particle_distributions[scalar_index_1(x, y, z, 1)] + particle_distributions[scalar_index_2(x, y, z, 2)] + particle_distributions[scalar_index_3(x, y, z, 3)] +
                             particle_distributions[scalar_index_4(x, y, z, 4)] + particle_distributions[scalar_index_5(x, y, z, 5)] + particle_distributions[scalar_index_6(x, y, z, 6)] + particle_distributions[scalar_index_7(x, y, z, 7)] +
                             particle_distributions[scalar_index_8(x, y, z, 8)] + particle_distributions[scalar_index_9(x, y, z, 9)] + particle_distributions[scalar_index_10(x, y, z, 10)] + particle_distributions[scalar_index_11(x, y, z, 11)] +
                             particle_distributions[scalar_index_12(x, y, z, 12)] + particle_distributions[scalar_index_13(x, y, z, 13)] + particle_distributions[scalar_index_14(x, y, z, 14)];

        double x_sum = particle_distributions[scalar_index_1(x, y, z, 1)] - particle_distributions[scalar_index_2(x, y, z, 2)] + particle_distributions[scalar_index_7(x, y, z, 7)] - particle_distributions[scalar_index_8(x, y, z, 8)] +
                       particle_distributions[scalar_index_9(x, y, z, 9)] - particle_distributions[scalar_index_10(x, y, z, 10)] + particle_distributions[scalar_index_11(x, y, z, 11)] - particle_distributions[scalar_index_12(x, y, z, 12)] -
                       particle_distributions[scalar_index_13(x, y, z, 13)] + particle_distributions[scalar_index_14(x, y, z, 14)];

        double y_sum = particle_distributions[scalar_index_3(x, y, z, 3)] - particle_distributions[scalar_index_4(x, y, z, 4)] + particle_distributions[scalar_index_7(x, y, z, 7)] - particle_distributions[scalar_index_8(x, y, z, 8)] +
                       particle_distributions[scalar_index_9(x, y, z, 9)] - particle_distributions[scalar_index_10(x, y, z, 10)] - particle_distributions[scalar_index_11(x, y, z, 11)] + particle_distributions[scalar_index_12(x, y, z, 12)] +
                       particle_distributions[scalar_index_13(x, y, z, 13)] - particle_distributions[scalar_index_14(x, y, z, 14)];

        double z_sum = particle_distributions[scalar_index_5(x, y, z, 5)] - particle_distributions[scalar_index_6(x, y, z, 6)] + particle_distributions[scalar_index_7(x, y, z, 7)] - particle_distributions[scalar_index_8(x, y, z, 8)] -
                       particle_distributions[scalar_index_9(x, y, z, 9)] + particle_distributions[scalar_index_10(x, y, z, 10)] + particle_distributions[scalar_index_11(x, y, z, 11)] - particle_distributions[scalar_index_12(x, y, z, 12)] +
                       particle_distributions[scalar_index_13(x, y, z, 13)] - particle_distributions[scalar_index_14(x, y, z, 14)];

        density_field[scalar_index(x, y, z)] = new_density;
        velocity_field_x[scalar_index(x, y, z)] = x_sum / new_density;
        velocity_field_y[scalar_index(x, y, z)] = y_sum / new_density;
        velocity_field_z[scalar_index(x, y, z)] = z_sum / new_density;
      }
    }
  }
}

void stream() {}

#define c_s_2_inv (1 / c_s_2)
#define c_s_4_inv (1 / c_s_4)
#define tau_weights_29 (weights_29 * tauinv)
#define tau_weights_19 (weights_19 * tauinv)
#define tau_weights_172 (weights_172 * tauinv)
#define c_s_2_2_inv (1 / c_s_2_2)
#define tauinv (1.0 / tau)
#define omtauinv (1.0 - tauinv)
void collision() { // Performs the collision step.
  for (int z = 0; z < NZ; z++) {
    for (int y = 0; y < NY; y++) {
      for (int x = 0; x < NX; x++) {
        double norm_square = ((velocity_field_x[scalar_index(x, y, z)]) * (velocity_field_x[scalar_index(x, y, z)]) + velocity_field_y[scalar_index(x, y, z)] * velocity_field_y[scalar_index(x, y, z)] +
                              velocity_field_z[scalar_index(x, y, z)] * velocity_field_z[scalar_index(x, y, z)]) *
                             c_s_2_2_inv;
        double norm_comp = 1.0 - norm_square;
        double curr_dens_field = density_field[scalar_index(x, y, z)];

        double feq = tau_weights_29 * curr_dens_field * norm_comp;
        particle_distributions[scalar_index(x, y, z, 0)] = omtauinv * particle_distributions[scalar_index(x, y, z, 0)] + feq;

        double feq1 = tau_weights_19 * curr_dens_field * (velocity_field_x[scalar_index(x, y, z)] * c_s_2_inv + velocity_field_x[scalar_index(x, y, z)] * velocity_field_x[scalar_index(x, y, z)] * c_s_4_inv + norm_comp);
        particle_distributions[scalar_index_1(x, y, z, 1)] = omtauinv * particle_distributions[scalar_index_1(x, y, z, 1)] + feq1;

        double feq2 = tau_weights_19 * curr_dens_field * (-velocity_field_x[scalar_index(x, y, z)] * c_s_2_inv + velocity_field_x[scalar_index(x, y, z)] * velocity_field_x[scalar_index(x, y, z)] * c_s_4_inv + norm_comp);
        particle_distributions[scalar_index_2(x, y, z, 2)] = omtauinv * particle_distributions[scalar_index_2(x, y, z, 2)] + feq2;

        double feq3 = tau_weights_19 * curr_dens_field * (velocity_field_y[scalar_index(x, y, z)] * c_s_2_inv + velocity_field_y[scalar_index(x, y, z)] * velocity_field_y[scalar_index(x, y, z)] * c_s_4_inv + norm_comp);
        particle_distributions[scalar_index_3(x, y, z, 3)] = omtauinv * particle_distributions[scalar_index_3(x, y, z, 3)] + feq3;

        double feq4 = tau_weights_19 * curr_dens_field * (-velocity_field_y[scalar_index(x, y, z)] * c_s_2_inv + velocity_field_y[scalar_index(x, y, z)] * velocity_field_y[scalar_index(x, y, z)] * c_s_4_inv + norm_comp);
        particle_distributions[scalar_index_4(x, y, z, 4)] = omtauinv * particle_distributions[scalar_index_4(x, y, z, 4)] + feq4;

        double feq5 = tau_weights_19 * curr_dens_field * (velocity_field_z[scalar_index(x, y, z)] * c_s_2_inv + velocity_field_z[scalar_index(x, y, z)] * velocity_field_z[scalar_index(x, y, z)] * c_s_4_inv + norm_comp);
        particle_distributions[scalar_index_5(x, y, z, 5)] = omtauinv * particle_distributions[scalar_index_5(x, y, z, 5)] + feq5;

        double feq6 = tau_weights_19 * curr_dens_field * (-velocity_field_z[scalar_index(x, y, z)] * c_s_2_inv + velocity_field_z[scalar_index(x, y, z)] * velocity_field_z[scalar_index(x, y, z)] * c_s_4_inv + norm_comp);
        particle_distributions[scalar_index_6(x, y, z, 6)] = omtauinv * particle_distributions[scalar_index_6(x, y, z, 6)] + feq6;

        double dot_product7 = velocity_field_x[scalar_index(x, y, z)] + velocity_field_y[scalar_index(x, y, z)] + velocity_field_z[scalar_index(x, y, z)];
        double feq7 = tau_weights_172 * curr_dens_field * (dot_product7 * c_s_2_inv + dot_product7 * dot_product7 * c_s_4_inv + norm_comp);
        particle_distributions[scalar_index_7(x, y, z, 7)] = omtauinv * particle_distributions[scalar_index_7(x, y, z, 7)] + feq7;

        double dot_product8 = -velocity_field_x[scalar_index(x, y, z)] - velocity_field_y[scalar_index(x, y, z)] - velocity_field_z[scalar_index(x, y, z)];
        double feq8 = tau_weights_172 * curr_dens_field * (dot_product8 * c_s_2_inv + dot_product8 * dot_product8 * c_s_4_inv + norm_comp);
        particle_distributions[scalar_index_8(x, y, z, 8)] = omtauinv * particle_distributions[scalar_index_8(x, y, z, 8)] + feq8;

        double dot_product9 = velocity_field_x[scalar_index(x, y, z)] + velocity_field_y[scalar_index(x, y, z)] - velocity_field_z[scalar_index(x, y, z)];
        double feq9 = tau_weights_172 * curr_dens_field * (dot_product9 * c_s_2_inv + dot_product9 * dot_product9 * c_s_4_inv + norm_comp);
        particle_distributions[scalar_index_9(x, y, z, 9)] = omtauinv * particle_distributions[scalar_index_9(x, y, z, 9)] + feq9;

        double dot_product10 = -velocity_field_x[scalar_index(x, y, z)] - velocity_field_y[scalar_index(x, y, z)] + velocity_field_z[scalar_index(x, y, z)];
        double feq10 = tau_weights_172 * curr_dens_field * (dot_product10 * c_s_2_inv + dot_product10 * dot_product10 * c_s_4_inv + norm_comp);
        particle_distributions[scalar_index_10(x, y, z, 10)] = omtauinv * particle_distributions[scalar_index_10(x, y, z, 10)] + feq10;

        double dot_product11 = velocity_field_x[scalar_index(x, y, z)] - velocity_field_y[scalar_index(x, y, z)] + velocity_field_z[scalar_index(x, y, z)];
        double feq11 = tau_weights_172 * curr_dens_field * (dot_product11 * c_s_2_inv + dot_product11 * dot_product11 * c_s_4_inv + norm_comp);
        particle_distributions[scalar_index_11(x, y, z, 11)] = omtauinv * particle_distributions[scalar_index_11(x, y, z, 11)] + feq11;

        double dot_product12 = -velocity_field_x[scalar_index(x, y, z)] + velocity_field_y[scalar_index(x, y, z)] - velocity_field_z[scalar_index(x, y, z)];
        double feq12 = tau_weights_172 * curr_dens_field * (dot_product12 * c_s_2_inv + dot_product12 * dot_product12 * c_s_4_inv + norm_comp);
        particle_distributions[scalar_index_12(x, y, z, 12)] = omtauinv * particle_distributions[scalar_index_12(x, y, z, 12)] + feq12;

        double dot_product13 = -velocity_field_x[scalar_index(x, y, z)] + velocity_field_y[scalar_index(x, y, z)] + velocity_field_z[scalar_index(x, y, z)];
        double feq13 = tau_weights_172 * curr_dens_field * (dot_product13 * c_s_2_inv + dot_product13 * dot_product13 * c_s_4_inv + norm_comp);
        particle_distributions[scalar_index_13(x, y, z, 13)] = omtauinv * particle_distributions[scalar_index_13(x, y, z, 13)] + feq13;

        double dot_product14 = velocity_field_x[scalar_index(x, y, z)] - velocity_field_y[scalar_index(x, y, z)] - velocity_field_z[scalar_index(x, y, z)];
        double feq14 = tau_weights_172 * curr_dens_field * (dot_product14 * c_s_2_inv + dot_product14 * dot_product14 * c_s_4_inv + norm_comp);
        particle_distributions[scalar_index_14(x, y, z, 14)] = omtauinv * particle_distributions[scalar_index_14(x, y, z, 14)] + feq14;
      }
    }
  }
}

void perform_timestep() {

// ----------------- COLLISION -----------------
#ifdef BENCHMARK
  if (PAPI_start(papi_event_set) != PAPI_OK) {
    fprintf(stderr, "PAPI start error!\n");
    exit(1);
  }
#endif
  start_run(collision_profiler);
  collision();
  end_run(collision_profiler);
#ifdef BENCHMARK
  if (PAPI_stop(papi_event_set, papi_values) != PAPI_OK) {
    fprintf(stderr, "PAPI stop error!\n");
    exit(1);
  }
  papi_collision_values[0] += 64 * papi_values[0];
  papi_collision_values[1] += papi_values[1];
#endif
  time_lbm++;
  time_lbm_x = time_lbm % NX;
  time_lbm_y = time_lbm % NY;
  time_lbm_z = time_lbm % NZ;
// // ----------------- STREAM -----------------
// #ifdef BENCHMARK
//   if (PAPI_start(papi_event_set) != PAPI_OK) {
//     fprintf(stderr, "PAPI start error!\n");
//     exit(1);
//   }
// #endif
//   start_run(stream_profiler);
//   stream();
//   end_run(stream_profiler);
// #ifdef BENCHMARK
//   if (PAPI_stop(papi_event_set, papi_values) != PAPI_OK) {
//     fprintf(stderr, "PAPI stop error!\n");
//     exit(1);
//   }
//   papi_stream_values[0] += 64 * papi_values[0];
//   papi_stream_values[1] += papi_values[1];
// #endif
// ----------------- DENSITY_MOMENTUM_MOMENT -----------------
#ifdef BENCHMARK
  if (PAPI_start(papi_event_set) != PAPI_OK) {
    fprintf(stderr, "PAPI start error!\n");
    exit(1);
  }
#endif
  start_run(compute_density_momentum_profiler);
  compute_density_momentum_moment();
  end_run(compute_density_momentum_profiler);
#ifdef BENCHMARK
  if (PAPI_stop(papi_event_set, papi_values) != PAPI_OK) {
    fprintf(stderr, "PAPI stop error!\n");
    exit(1);
  }
  papi_density_values[0] += 64 * papi_values[0];
  papi_density_values[1] += papi_values[1];
#endif
}

int main(int argc, char const *argv[]) {
  if (system("rm -rf output") == -1) {
    std::cout << "UNSUCCESSFULLY REMOVED OUTPUT FOLDER" << std::endl;
  }
  if (system("mkdir output") == -1) {
    std::cout << "UNSUCCESSFULLY CREATED OUTPUT FOLDER" << std::endl;
  }

  unsigned long long start_cycle, end_cycle;
  time_t start_sec, end_sec;

#ifdef TIMING
  asm volatile("RDTSC" : "=A"(start_cycle));
  time(&start_sec);
#endif

  // setup LBM
  initialise();

  double viscosity = c_s * c_s * (tau - 0.5);

#ifdef OUTPUT
  output_lbm_data("output/0.csv", true);
  output_indices_file();
#endif

  int runs = NT;

#ifdef PROFILE
  const int PROFILE_RUNS = 5;
  const int PROFILE_DIGITS = floor(log10(PROFILE_RUNS)) + 1;
  printf("\rRun %-*d/%d done", PROFILE_DIGITS, 0, PROFILE_RUNS);
  fflush(stdout);
  for (int profile_counter = 0; profile_counter < PROFILE_RUNS; profile_counter++) {
    time_lbm = 0;
#endif

    // finish setup

    // start simulation
    for (int i = 0; i < runs; i = i + 1) {
      perform_timestep();

#ifndef TIMING
#ifdef OUTPUT
      if ((i + 1) % save_every == 0) {
        double percentage = (double)(i + 1) / (double)(runs)*100.0;
        std::cout << "Saving data - " << (i + 1) << "/" << runs << " (" << percentage << "%)" << '\n';
        output_lbm_data("output/" + std::to_string(i + 1) + ".csv", true);
      }
      std::ofstream timestamp_file;
      timestamp_file.open("output/timestamp_" + std::to_string(NX) + "_" + std::to_string(NY) + "_" + std::to_string(NZ) + "_" + std::to_string(NT) + ".txt");
      timestamp_file.close();
#endif
#endif
    }
    std::cout << std::endl;
#ifdef PROFILE
    printf("\rRun %-*d/%d done", PROFILE_DIGITS, profile_counter + 1, PROFILE_RUNS);
    fflush(stdout);
  }
#endif
  free_up();

#ifdef TIMING
  asm volatile("RDTSC" : "=A"(end_cycle));
  time(&end_sec);
  printf("Cycles taken: %llu (%ld seconds)\n", end_cycle - start_cycle, end_sec - start_sec);

#endif

#ifdef PROFILE
  printf("\nProfiling results:\n");
  profiler_stats compute_density_stats = finish_profiler(compute_density_momentum_profiler);
  printf("- Compute density  Calculation: %4.2f Flops/Cycle, %10ld cycles in %d runs. "
         "Arithmetic intensity: %4.2f\n",
         compute_density_stats.performance, compute_density_stats.cycles, compute_density_stats.runs, compute_density_stats.arithmetic_intensity);

  profiler_stats collision_stats = finish_profiler(collision_profiler);
  printf("- Compute collision  Calculation: %4.2f Flops/Cycle, %10ld cycles in %d runs. "
         "Arithmetic intensity: %4.2f\n",
         collision_stats.performance, collision_stats.cycles, collision_stats.runs, collision_stats.arithmetic_intensity);
  profiler_stats stream_stats = finish_profiler(stream_profiler);
  printf("- Compute stream  Calculation: %4.2f Flops/Cycle, %10ld cycles in %d runs. "
         "Arithmetic intensity: %4.2f\n",
         stream_stats.performance, stream_stats.cycles, stream_stats.runs, stream_stats.arithmetic_intensity);
#endif

#ifdef BENCHMARK
  printf("Density Momentum Moment Mem Transfer: %lld\n", papi_density_values[0] / NT);
  printf("Density Momentum Moment Floating point operations: %lld\n", papi_density_values[1] / NT);
  printf("Density Momentum Moment Arithmetic Intensity: %f\n", (double)papi_density_values[1] / (double)papi_density_values[0]);

  printf("Collision Mem Transfer: %lld\n", papi_collision_values[0] / NT);
  printf("Collision Floating point operations: %lld\n", papi_collision_values[1] / NT);
  printf("Collision Arithmetic Intensity: %f\n", (double)papi_collision_values[1] / (double)papi_collision_values[0]);

  printf("Stream Mem Transfer: %lld\n", papi_stream_values[0] / NT);
  printf("Stream Floating point operations: %lld\n", papi_stream_values[1] / NT);
  printf("Stream Arithmetic Intensity: %f\n", (double)papi_stream_values[1] / (double)papi_stream_values[0]);

  // Cleanup PAPI
  PAPI_shutdown();
#endif

  return 0;
}
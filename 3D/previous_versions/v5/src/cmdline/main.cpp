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
double *velocity_field;
double *previous_particle_distributions;
double *particle_distributions;
int direction_size = 15;

#ifdef BENCHMARK
int papi_event_set = PAPI_NULL;

long long papi_values[2];

long long papi_density_values[2] = {0, 0};
long long papi_collision_values[2] = {0, 0};
long long papi_stream_values[2] = {0, 0};
#endif

inline int scalar_index(int x, int y, int z) { return (z * NX * NY) + (y * NX) + x; }
inline int scalar_index(int x, int y, int z, int w) { return (x + y * NX + z * NX * NY + w * NX * NY * NZ); }

// Lattice directions using D3DQ15. assumed speed of sound c_s = 1/sqrt(3).
vector_3_int *directions;
double *weights;
int *reverse_indexes;

void initialise() {
  int box_flatten_length = NX * NY * NZ;
  int distributions_flatten_length = box_flatten_length * direction_size;
  density_field = (double *)malloc(box_flatten_length * sizeof(double));
  velocity_field = (double *)malloc(box_flatten_length * 3 * sizeof(double));

  previous_particle_distributions = (double *)malloc(distributions_flatten_length * sizeof(double));
  particle_distributions = (double *)malloc(distributions_flatten_length * sizeof(double));
  srand(42);

  for (int i = 0; i < NX * NY * NZ; i++) {
    density_field[i] = 1.0;
  }
  for (int x = 0; x < NX; x++) {
    for (int y = 0; y < NY; y++) {
      for (int z = 0; z < NZ; z++) {
        for (int i = 0; i < direction_size; i++) {
          previous_particle_distributions[scalar_index(x, y, z, i)] = weights[i];
          particle_distributions[scalar_index(x, y, z, i)] = weights[i];
        }
        velocity_field[scalar_index(x, y, z) * 3] = (double)(rand());
        velocity_field[scalar_index(x, y, z) * 3 + 1] = (double)(rand());
        velocity_field[scalar_index(x, y, z) * 3 + 2] = (double)(rand());
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

#define c_s_4 (2 * c_s * c_s * c_s * c_s)
#define c_s_2 (c_s * c_s)
#define c_s_2_2 (2 * c_s * c_s)
#define weights_172 (1.0 / 72.0)
#define weights_19 (1.0 / 9.0)
#define weights_19 (1.0 / 9.0)
#define weights_29 (2.0 / 9.0)

void set_velocity_set() {
  // Allocate memory for an array of vector_3_int structs
  directions = (vector_3_int *)malloc(15 * sizeof(vector_3_int));

  // Initialize each element of the array
  directions[0] = {0, 0, 0};
  directions[1] = {1, 0, 0};
  directions[2] = {-1, 0, 0};
  directions[3] = {0, 1, 0};
  directions[4] = {0, -1, 0};

  directions[5] = {0, 0, 1};
  directions[6] = {0, 0, -1};
  directions[7] = {1, 1, 1};
  directions[8] = {-1, -1, -1};
  directions[9] = {1, 1, -1};
  directions[10] = {-1, -1, 1};
  directions[11] = {1, -1, 1};
  directions[12] = {-1, 1, -1};
  directions[13] = {-1, 1, 1};
  directions[14] = {1, -1, -1};

  weights = (double *)malloc(15 * sizeof(double));
  // Initialize each element of the array
  weights[0] = 2.0 / 9.0;
  weights[1] = 1.0 / 9.0;
  weights[2] = 1.0 / 9.0;
  weights[3] = 1.0 / 9.0;
  weights[4] = 1.0 / 9.0;
  weights[5] = 1.0 / 9.0;
  weights[6] = 1.0 / 9.0;
  weights[7] = 1.0 / 72.0;
  weights[8] = 1.0 / 72.0;
  weights[9] = 1.0 / 72.0;
  weights[10] = 1.0 / 72.0;
  weights[11] = 1.0 / 72.0;
  weights[12] = 1.0 / 72.0;
  weights[13] = 1.0 / 72.0;
  weights[14] = 1.0 / 72.0;
}

void compute_density_momentum_moment() {
  int scalar_index_curr = 0;
  for (int z = 0; z < NZ; z++) {
    for (int y = 0; y < NY; y++) {
      for (int x = 0; x < NX; x++) {

        double new_density = 0, new_density_1 = 0, new_density_2 = 0, new_density_3 = 0, new_density_4 = 0, new_density_5 = 0;
        double x_1 = 0, x_2 = 0, x_3 = 0, x_4 = 0, x_5 = 0;
        double y_1 = 0, y_2 = 0, y_3 = 0, y_4 = 0, y_5 = 0;
        double z_1 = 0, z_2 = 0, z_3 = 0, z_4 = 0, z_5 = 0;
        double x_sum = 0.0;
        double y_sum = 0.0;
        double z_sum = 0.0;
        new_density += particle_distributions[scalar_index(x, y, z, 0)];

        new_density += particle_distributions[scalar_index(x, y, z, 1)];
        x_sum += particle_distributions[scalar_index(x, y, z, 1)];

        new_density += particle_distributions[scalar_index(x, y, z, 2)];
        x_sum -= particle_distributions[scalar_index(x, y, z, 2)];

        new_density += particle_distributions[scalar_index(x, y, z, 3)];
        y_sum += particle_distributions[scalar_index(x, y, z, 3)] * 1;

        new_density += particle_distributions[scalar_index(x, y, z, 4)];
        y_sum -= particle_distributions[scalar_index(x, y, z, 4)];

        new_density += particle_distributions[scalar_index(x, y, z, 5)];
        z_sum += particle_distributions[scalar_index(x, y, z, 5)];

        new_density += particle_distributions[scalar_index(x, y, z, 6)];
        z_sum -= particle_distributions[scalar_index(x, y, z, 6)];

        new_density += particle_distributions[scalar_index(x, y, z, 7)];
        x_sum += particle_distributions[scalar_index(x, y, z, 7)];
        y_sum += particle_distributions[scalar_index(x, y, z, 7)];
        z_sum += particle_distributions[scalar_index(x, y, z, 7)];

        new_density += particle_distributions[scalar_index(x, y, z, 8)];
        x_sum -= particle_distributions[scalar_index(x, y, z, 8)];
        y_sum -= particle_distributions[scalar_index(x, y, z, 8)];
        z_sum -= particle_distributions[scalar_index(x, y, z, 8)];

        new_density += particle_distributions[scalar_index(x, y, z, 9)];
        x_sum += particle_distributions[scalar_index(x, y, z, 9)];
        y_sum += particle_distributions[scalar_index(x, y, z, 9)];
        z_sum -= particle_distributions[scalar_index(x, y, z, 9)];

        new_density += particle_distributions[scalar_index(x, y, z, 10)];
        x_sum -= particle_distributions[scalar_index(x, y, z, 10)];
        y_sum -= particle_distributions[scalar_index(x, y, z, 10)];
        z_sum += particle_distributions[scalar_index(x, y, z, 10)];

        new_density += particle_distributions[scalar_index(x, y, z, 11)];
        x_sum += particle_distributions[scalar_index(x, y, z, 11)];
        y_sum -= particle_distributions[scalar_index(x, y, z, 11)];
        z_sum += particle_distributions[scalar_index(x, y, z, 11)];

        new_density += particle_distributions[scalar_index(x, y, z, 12)];
        x_sum -= particle_distributions[scalar_index(x, y, z, 12)];
        y_sum += particle_distributions[scalar_index(x, y, z, 12)];
        z_sum -= particle_distributions[scalar_index(x, y, z, 12)];

        new_density += particle_distributions[scalar_index(x, y, z, 13)];
        x_sum -= particle_distributions[scalar_index(x, y, z, 13)];
        y_sum += particle_distributions[scalar_index(x, y, z, 13)];
        z_sum += particle_distributions[scalar_index(x, y, z, 13)];

        new_density += particle_distributions[scalar_index(x, y, z, 14)];
        x_sum += particle_distributions[scalar_index(x, y, z, 14)];
        y_sum -= particle_distributions[scalar_index(x, y, z, 14)];
        z_sum -= particle_distributions[scalar_index(x, y, z, 14)];

        density_field[scalar_index(x, y, z)] = new_density;
        velocity_field[scalar_index(x, y, z) * 3] = x_sum / new_density;
        velocity_field[scalar_index(x, y, z) * 3 + 1] = y_sum / new_density;
        velocity_field[scalar_index(x, y, z) * 3 + 2] = z_sum / new_density;
      }
    }
  }
}

void stream() {
  for (int z = 0; z < NZ; z++) {
    for (int y = 0; y < NY; y++) {
      for (int x = 0; x < NX; x++) {
        int xmd = x;
        int ymd = y;
        int zmd = z;
        particle_distributions[scalar_index(x, y, z, 0)] = previous_particle_distributions[scalar_index(xmd, ymd, zmd, 0)];

        xmd = (NX + x - 1) % NX;
        ymd = y;
        zmd = z;
        particle_distributions[scalar_index(x, y, z, 1)] = previous_particle_distributions[scalar_index(xmd, ymd, zmd, 1)];

        xmd = (NX + x + 1) % NX;
        ymd = y;
        zmd = z;
        particle_distributions[scalar_index(x, y, z, 2)] = previous_particle_distributions[scalar_index(xmd, ymd, zmd, 2)];

        xmd = x;
        ymd = (NY + y - 1) % NY;
        zmd = z;
        particle_distributions[scalar_index(x, y, z, 3)] = previous_particle_distributions[scalar_index(xmd, ymd, zmd, 3)];

        xmd = x;
        ymd = (NY + y + 1) % NY;
        zmd = z;
        particle_distributions[scalar_index(x, y, z, 4)] = previous_particle_distributions[scalar_index(xmd, ymd, zmd, 4)];

        xmd = x;
        ymd = y;
        zmd = (NZ + z - 1) % NZ;
        particle_distributions[scalar_index(x, y, z, 5)] = previous_particle_distributions[scalar_index(xmd, ymd, zmd, 5)];

        xmd = x;
        ymd = y;
        zmd = (NZ + z + 1) % NZ;
        particle_distributions[scalar_index(x, y, z, 6)] = previous_particle_distributions[scalar_index(xmd, ymd, zmd, 6)];

        xmd = (NX + x - 1) % NX;
        ymd = (NY + y - 1) % NY;
        zmd = (NZ + z - 1) % NZ;
        particle_distributions[scalar_index(x, y, z, 7)] = previous_particle_distributions[scalar_index(xmd, ymd, zmd, 7)];

        xmd = (NX + x + 1) % NX;
        ymd = (NY + y + 1) % NY;
        zmd = (NZ + z + 1) % NZ;
        particle_distributions[scalar_index(x, y, z, 8)] = previous_particle_distributions[scalar_index(xmd, ymd, zmd, 8)];

        xmd = (NX + x - 1) % NX;
        ymd = (NY + y - 1) % NY;
        zmd = (NZ + z + 1) % NZ;
        particle_distributions[scalar_index(x, y, z, 9)] = previous_particle_distributions[scalar_index(xmd, ymd, zmd, 9)];

        xmd = (NX + x + 1) % NX;
        ymd = (NY + y + 1) % NY;
        zmd = (NZ + z - 1) % NZ;
        particle_distributions[scalar_index(x, y, z, 10)] = previous_particle_distributions[scalar_index(xmd, ymd, zmd, 10)];

        xmd = (NX + x - 1) % NX;
        ymd = (NY + y + 1) % NY;
        zmd = (NZ + z - 1) % NZ;
        particle_distributions[scalar_index(x, y, z, 11)] = previous_particle_distributions[scalar_index(xmd, ymd, zmd, 11)];

        xmd = (NX + x + 1) % NX;
        ymd = (NY + y - 1) % NY;
        zmd = (NZ + z + 1) % NZ;
        particle_distributions[scalar_index(x, y, z, 12)] = previous_particle_distributions[scalar_index(xmd, ymd, zmd, 12)];

        xmd = (NX + x + 1) % NX;
        ymd = (NY + y - 1) % NY;
        zmd = (NZ + z - 1) % NZ;
        particle_distributions[scalar_index(x, y, z, 13)] = previous_particle_distributions[scalar_index(xmd, ymd, zmd, 13)];

        xmd = (NX + x - 1) % NX;
        ymd = (NY + y + 1) % NY;
        zmd = (NZ + z + 1) % NZ;
        particle_distributions[scalar_index(x, y, z, 14)] = previous_particle_distributions[scalar_index(xmd, ymd, zmd, 14)];
      }
    }
  }
}
void collision() { // Performs the collision step.
  const double tauinv = 1.0 / tau;
  const double omtauinv = 1.0 - tauinv; // 1 - 1/tau
  for (int z = 0; z < NZ; z++) {
    for (int y = 0; y < NY; y++) {
      for (int x = 0; x < NX; x++) {
        double norm_square = ((velocity_field[scalar_index(x, y, z) * 3]) * (velocity_field[scalar_index(x, y, z) * 3]) + velocity_field[scalar_index(x, y, z) * 3 + 1] * velocity_field[scalar_index(x, y, z) * 3 + 1] +
                              velocity_field[scalar_index(x, y, z) * 3 + 2] * velocity_field[scalar_index(x, y, z) * 3 + 2]) /
                             c_s_2_2;

        double feq = weights_29 * density_field[scalar_index(x, y, z)] * (1.0 - norm_square);
        previous_particle_distributions[scalar_index(x, y, z, 0)] = omtauinv * particle_distributions[scalar_index(x, y, z, 0)] + tauinv * feq;

        double dot_product = velocity_field[scalar_index(x, y, z) * 3];
        feq = weights_19 * density_field[scalar_index(x, y, z)] * (1.0 + dot_product / c_s_2 + dot_product * dot_product / c_s_4 - norm_square);
        previous_particle_distributions[scalar_index(x, y, z, 1)] = omtauinv * particle_distributions[scalar_index(x, y, z, 1)] + tauinv * feq;

        dot_product = -velocity_field[scalar_index(x, y, z) * 3];
        feq = weights_19 * density_field[scalar_index(x, y, z)] * (1.0 + dot_product / c_s_2 + dot_product * dot_product / c_s_4 - norm_square);
        previous_particle_distributions[scalar_index(x, y, z, 2)] = omtauinv * particle_distributions[scalar_index(x, y, z, 2)] + tauinv * feq;

        dot_product = velocity_field[scalar_index(x, y, z) * 3 + 1] + velocity_field[scalar_index(x, y, z) * 3 + 2] * (double)0;
        feq = weights_19 * density_field[scalar_index(x, y, z)] * (1.0 + dot_product / c_s_2 + dot_product * dot_product / c_s_4 - norm_square);
        previous_particle_distributions[scalar_index(x, y, z, 3)] = omtauinv * particle_distributions[scalar_index(x, y, z, 3)] + tauinv * feq;

        dot_product = -velocity_field[scalar_index(x, y, z) * 3 + 1];
        feq = weights_19 * density_field[scalar_index(x, y, z)] * (1.0 + dot_product / c_s_2 + dot_product * dot_product / c_s_4 - norm_square);
        previous_particle_distributions[scalar_index(x, y, z, 4)] = omtauinv * particle_distributions[scalar_index(x, y, z, 4)] + tauinv * feq;

        dot_product = velocity_field[scalar_index(x, y, z) * 3] * (double)0 + velocity_field[scalar_index(x, y, z) * 3 + 1] * (double)0 + velocity_field[scalar_index(x, y, z) * 3 + 2] * (double)1;
        feq = weights_19 * density_field[scalar_index(x, y, z)] * (1.0 + dot_product / c_s_2 + dot_product * dot_product / c_s_4 - norm_square);
        previous_particle_distributions[scalar_index(x, y, z, 5)] = omtauinv * particle_distributions[scalar_index(x, y, z, 5)] + tauinv * feq;

        dot_product = velocity_field[scalar_index(x, y, z) * 3] * (double)0 + velocity_field[scalar_index(x, y, z) * 3 + 1] * (double)0 + velocity_field[scalar_index(x, y, z) * 3 + 2] * (double)-1;
        feq = weights_19 * density_field[scalar_index(x, y, z)] * (1.0 + dot_product / c_s_2 + dot_product * dot_product / c_s_4 - norm_square);
        previous_particle_distributions[scalar_index(x, y, z, 6)] = omtauinv * particle_distributions[scalar_index(x, y, z, 6)] + tauinv * feq;

        dot_product = velocity_field[scalar_index(x, y, z) * 3] * (double)1 + velocity_field[scalar_index(x, y, z) * 3 + 1] * (double)1 + velocity_field[scalar_index(x, y, z) * 3 + 2] * (double)1;
        feq = weights_172 * density_field[scalar_index(x, y, z)] * (1.0 + dot_product / c_s_2 + dot_product * dot_product / c_s_4 - norm_square);
        previous_particle_distributions[scalar_index(x, y, z, 7)] = omtauinv * particle_distributions[scalar_index(x, y, z, 7)] + tauinv * feq;

        dot_product = velocity_field[scalar_index(x, y, z) * 3] * (double)-1 + velocity_field[scalar_index(x, y, z) * 3 + 1] * (double)-1 + velocity_field[scalar_index(x, y, z) * 3 + 2] * (double)-1;
        feq = weights_172 * density_field[scalar_index(x, y, z)] * (1.0 + dot_product / c_s_2 + dot_product * dot_product / c_s_4 - norm_square);
        previous_particle_distributions[scalar_index(x, y, z, 8)] = omtauinv * particle_distributions[scalar_index(x, y, z, 8)] + tauinv * feq;

        dot_product = velocity_field[scalar_index(x, y, z) * 3] * (double)1 + velocity_field[scalar_index(x, y, z) * 3 + 1] * (double)1 + velocity_field[scalar_index(x, y, z) * 3 + 2] * (double)-1;
        feq = weights_172 * density_field[scalar_index(x, y, z)] * (1.0 + dot_product / c_s_2 + dot_product * dot_product / c_s_4 - norm_square);
        previous_particle_distributions[scalar_index(x, y, z, 9)] = omtauinv * particle_distributions[scalar_index(x, y, z, 9)] + tauinv * feq;

        dot_product = -velocity_field[scalar_index(x, y, z) * 3] - velocity_field[scalar_index(x, y, z) * 3 + 1] + velocity_field[scalar_index(x, y, z) * 3 + 2] * (double)1;
        feq = weights_172 * density_field[scalar_index(x, y, z)] * (1.0 + dot_product / c_s_2 + dot_product * dot_product / c_s_4 - norm_square);
        previous_particle_distributions[scalar_index(x, y, z, 10)] = omtauinv * particle_distributions[scalar_index(x, y, z, 10)] + tauinv * feq;

        dot_product = velocity_field[scalar_index(x, y, z) * 3] - velocity_field[scalar_index(x, y, z) * 3 + 1] + velocity_field[scalar_index(x, y, z) * 3 + 2];
        feq = weights_172 * density_field[scalar_index(x, y, z)] * (1.0 + dot_product / c_s_2 + dot_product * dot_product / c_s_4 - norm_square);
        previous_particle_distributions[scalar_index(x, y, z, 11)] = omtauinv * particle_distributions[scalar_index(x, y, z, 11)] + tauinv * feq;

        dot_product = -velocity_field[scalar_index(x, y, z) * 3] + velocity_field[scalar_index(x, y, z) * 3 + 1] - velocity_field[scalar_index(x, y, z) * 3 + 2];
        feq = weights_172 * density_field[scalar_index(x, y, z)] * (1.0 + dot_product / c_s_2 + dot_product * dot_product / c_s_4 - norm_square);
        previous_particle_distributions[scalar_index(x, y, z, 12)] = omtauinv * particle_distributions[scalar_index(x, y, z, 12)] + tauinv * feq;

        dot_product = -velocity_field[scalar_index(x, y, z) * 3] + velocity_field[scalar_index(x, y, z) * 3 + 1] + velocity_field[scalar_index(x, y, z) * 3 + 2] * (double)1;
        feq = weights_172 * density_field[scalar_index(x, y, z)] * (1.0 + dot_product / c_s_2 + dot_product * dot_product / c_s_4 - norm_square);
        previous_particle_distributions[scalar_index(x, y, z, 13)] = omtauinv * particle_distributions[scalar_index(x, y, z, 13)] + tauinv * feq;

        dot_product = velocity_field[scalar_index(x, y, z) * 3] - velocity_field[scalar_index(x, y, z) * 3 + 1] - velocity_field[scalar_index(x, y, z) * 3 + 2];
        feq = weights_172 * density_field[scalar_index(x, y, z)] * (1.0 + dot_product / c_s_2 + dot_product * dot_product / c_s_4 - norm_square);
        previous_particle_distributions[scalar_index(x, y, z, 14)] = omtauinv * particle_distributions[scalar_index(x, y, z, 14)] + tauinv * feq;
      }
    }
  }
}

void perform_timestep() {
  time_lbm++;

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

// ----------------- STREAM -----------------
#ifdef BENCHMARK
  if (PAPI_start(papi_event_set) != PAPI_OK) {
    fprintf(stderr, "PAPI start error!\n");
    exit(1);
  }
#endif
  start_run(stream_profiler);
  stream();
  end_run(stream_profiler);
#ifdef BENCHMARK
  if (PAPI_stop(papi_event_set, papi_values) != PAPI_OK) {
    fprintf(stderr, "PAPI stop error!\n");
    exit(1);
  }
  papi_stream_values[0] += 64 * papi_values[0];
  papi_stream_values[1] += papi_values[1];
#endif
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
  set_velocity_set();
  initialise();

  double viscosity = c_s * c_s * (tau - 0.5);
  // std::cout<<viscosity<<std::endl;

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
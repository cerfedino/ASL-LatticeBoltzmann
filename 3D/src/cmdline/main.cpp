#include <iostream>
#include <stdio.h>
// Now Linux only.
#include "profile.h"
#include "utils.h"
#include <fstream>
#include <immintrin.h>
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
double *velocity_field_x;
double *velocity_field_y;
double *velocity_field_z;
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
  density_field = (double *)aligned_alloc(32, box_flatten_length * sizeof(double));
  velocity_field_x = (double *)aligned_alloc(32, box_flatten_length * sizeof(double));
  velocity_field_y = (double *)aligned_alloc(32, box_flatten_length * sizeof(double));
  velocity_field_z = (double *)aligned_alloc(32, box_flatten_length * sizeof(double));
  previous_particle_distributions = (double *)aligned_alloc(32, distributions_flatten_length * sizeof(double));
  particle_distributions = (double *)aligned_alloc(32, distributions_flatten_length * sizeof(double));
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

        velocity_field_x[scalar_index(x, y, z)] = 0.5;
        velocity_field_y[scalar_index(x, y, z)] = 0.25;
        velocity_field_z[scalar_index(x, y, z)] = 0.25;
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
#define c_s_4_inv (1.0 / (2 * c_s * c_s * c_s * c_s))
#define c_s_2 (c_s * c_s)
#define c_s_2_inv (1.0 / (c_s * c_s))
#define min_c_s_2_inv (-1.0 / (c_s * c_s))
#define c_s_2_2 (2 * c_s * c_s)
#define c_s_2_2_inv (1.0 / (2 * c_s * c_s))
#define weights_172 (1.0 / 72.0)
#define weights_19 (1.0 / 9.0)
#define weights_19 (1.0 / 9.0)
#define weights_29 (2.0 / 9.0)
#define tauinv (1.0 / tau)
#define omtauinv (1.0 - tauinv)

__m256d const_omtauinv_vec = _mm256_set1_pd(omtauinv);
__m256d const_tauinv_vec = _mm256_set1_pd(tauinv);
__m256d const_cs_2_2_vec = _mm256_set1_pd(c_s_2_2);
__m256d const_cs_2_2_inv_vec = _mm256_set1_pd(c_s_2_2_inv);
__m256d const_cs_2_inv_vec = _mm256_set1_pd(c_s_2_inv);
__m256d const_min_cs_2_inv_vec = _mm256_set1_pd(min_c_s_2_inv);

__m256d const_cs_4_inv_vec = _mm256_set1_pd(c_s_4_inv);
__m256d const_1_vec = _mm256_set1_pd(1.0);
__m256d const_weights_29_vec = _mm256_set1_pd(weights_29);
__m256d const_weights_19_vec = _mm256_set1_pd(weights_19);
__m256d const_weights_172_vec = _mm256_set1_pd(weights_172);

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

  for (int z = 0; z < NZ; z++) {
    for (int y = 0; y < NY; y++) {
      for (int x = 0; x < NX; x++) {
        double new_density = particle_distributions[scalar_index(x, y, z, 0)] + particle_distributions[scalar_index(x, y, z, 1)] + particle_distributions[scalar_index(x, y, z, 2)] + particle_distributions[scalar_index(x, y, z, 3)] +
                             particle_distributions[scalar_index(x, y, z, 4)] + particle_distributions[scalar_index(x, y, z, 5)] + particle_distributions[scalar_index(x, y, z, 6)] + particle_distributions[scalar_index(x, y, z, 7)] +
                             particle_distributions[scalar_index(x, y, z, 8)] + particle_distributions[scalar_index(x, y, z, 9)] + particle_distributions[scalar_index(x, y, z, 10)] + particle_distributions[scalar_index(x, y, z, 11)] +
                             particle_distributions[scalar_index(x, y, z, 12)] + particle_distributions[scalar_index(x, y, z, 13)] + particle_distributions[scalar_index(x, y, z, 14)];

        double x_sum = particle_distributions[scalar_index(x, y, z, 1)] - particle_distributions[scalar_index(x, y, z, 2)] + particle_distributions[scalar_index(x, y, z, 7)] - particle_distributions[scalar_index(x, y, z, 8)] +
                       particle_distributions[scalar_index(x, y, z, 9)] - particle_distributions[scalar_index(x, y, z, 10)] + particle_distributions[scalar_index(x, y, z, 11)] - particle_distributions[scalar_index(x, y, z, 12)] -
                       particle_distributions[scalar_index(x, y, z, 13)] + particle_distributions[scalar_index(x, y, z, 14)];

        double y_sum = particle_distributions[scalar_index(x, y, z, 3)] - particle_distributions[scalar_index(x, y, z, 4)] + particle_distributions[scalar_index(x, y, z, 7)] - particle_distributions[scalar_index(x, y, z, 8)] +
                       particle_distributions[scalar_index(x, y, z, 9)] - particle_distributions[scalar_index(x, y, z, 10)] - particle_distributions[scalar_index(x, y, z, 11)] + particle_distributions[scalar_index(x, y, z, 12)] +
                       particle_distributions[scalar_index(x, y, z, 13)] - particle_distributions[scalar_index(x, y, z, 14)];

        double z_sum = particle_distributions[scalar_index(x, y, z, 5)] - particle_distributions[scalar_index(x, y, z, 6)] + particle_distributions[scalar_index(x, y, z, 7)] - particle_distributions[scalar_index(x, y, z, 8)] -
                       particle_distributions[scalar_index(x, y, z, 9)] + particle_distributions[scalar_index(x, y, z, 10)] + particle_distributions[scalar_index(x, y, z, 11)] - particle_distributions[scalar_index(x, y, z, 12)] +
                       particle_distributions[scalar_index(x, y, z, 13)] - particle_distributions[scalar_index(x, y, z, 14)];

        density_field[scalar_index(x, y, z)] = new_density;
        velocity_field_x[scalar_index(x, y, z)] = x_sum / new_density;
        velocity_field_y[scalar_index(x, y, z)] = y_sum / new_density;
        velocity_field_z[scalar_index(x, y, z)] = z_sum / new_density;
      }
    }
  }
}

void stream() {
  int z = 0, y = 0, x = 0;

  for (z = 0; z < NZ; z++) {
    for (y = 0; y < NY; y++) {
      for (x = 0; x < NX; x++) {
        particle_distributions[scalar_index(x, y, z, 1)] = previous_particle_distributions[scalar_index((NX + x - 1) % NX, y, z, 1)];
        particle_distributions[scalar_index(x, y, z, 2)] = previous_particle_distributions[scalar_index((NX + x + 1) % NX, y, z, 2)];
        particle_distributions[scalar_index(x, y, z, 3)] = previous_particle_distributions[scalar_index(x, (NY + y - 1) % NY, z, 3)];
        particle_distributions[scalar_index(x, y, z, 4)] = previous_particle_distributions[scalar_index(x, (NY + y + 1) % NY, z, 4)];
        particle_distributions[scalar_index(x, y, z, 5)] = previous_particle_distributions[scalar_index(x, y, (NZ + z - 1) % NZ, 5)];
        particle_distributions[scalar_index(x, y, z, 6)] = previous_particle_distributions[scalar_index(x, y, (NZ + z + 1) % NZ, 6)];
        particle_distributions[scalar_index(x, y, z, 7)] = previous_particle_distributions[scalar_index((NX + x - 1) % NX, (NY + y - 1) % NY, (NZ + z - 1) % NZ, 7)];
        particle_distributions[scalar_index(x, y, z, 8)] = previous_particle_distributions[scalar_index((NX + x + 1) % NX, (NY + y + 1) % NY, (NZ + z + 1) % NZ, 8)];
        particle_distributions[scalar_index(x, y, z, 9)] = previous_particle_distributions[scalar_index((NX + x - 1) % NX, (NY + y - 1) % NY, (NZ + z + 1) % NZ, 9)];
        particle_distributions[scalar_index(x, y, z, 10)] = previous_particle_distributions[scalar_index((NX + x + 1) % NX, (NY + y + 1) % NY, (NZ + z - 1) % NZ, 10)];
        particle_distributions[scalar_index(x, y, z, 11)] = previous_particle_distributions[scalar_index((NX + x - 1) % NX, (NY + y + 1) % NY, (NZ + z - 1) % NZ, 11)];
        particle_distributions[scalar_index(x, y, z, 12)] = previous_particle_distributions[scalar_index((NX + x + 1) % NX, (NY + y - 1) % NY, (NZ + z + 1) % NZ, 12)];
        particle_distributions[scalar_index(x, y, z, 13)] = previous_particle_distributions[scalar_index((NX + x + 1) % NX, (NY + y - 1) % NY, (NZ + z - 1) % NZ, 13)];
        particle_distributions[scalar_index(x, y, z, 14)] = previous_particle_distributions[scalar_index((NX + x - 1) % NX, (NY + y + 1) % NY, (NZ + z + 1) % NZ, 14)];
      }
    }
  }
}

void collision() { // Performs the collision step.
  for (int z = 0; z < NZ; z++) {
    for (int y = 0; y < NY; y++) {
      for (int x = 0; x < NX; x += 4) {
        __m256d density_field_vec = _mm256_load_pd(&density_field[scalar_index(x, y, z)]);         // density_field[x, y, z]
        __m256d weights_19_mul_density = _mm256_mul_pd(const_weights_19_vec, density_field_vec);   // weights_19 * density_field[x, y, z]
        __m256d weights_172_mul_density = _mm256_mul_pd(const_weights_172_vec, density_field_vec); // weights_127 * density_field[x, y, z]

        __m256d velocity_x_vec = _mm256_load_pd(&velocity_field_x[scalar_index(x, y, z)]); // velocity_field_x[x, y, z]
        __m256d velocity_x_pow_vec = _mm256_mul_pd(velocity_x_vec, velocity_x_vec);        // velocity_field_x[x, y, z] ** 2
        __m256d velocity_y_vec = _mm256_load_pd(&velocity_field_y[scalar_index(x, y, z)]); // velocity_field_y[x, y, z]
        __m256d velocity_y_pow_vec = _mm256_mul_pd(velocity_y_vec, velocity_y_vec);        // velocity_field_y[x, y, z] ** 2
        __m256d velocity_z_vec = _mm256_load_pd(&velocity_field_z[scalar_index(x, y, z)]); // velocity_field_z[x, y, z]
        __m256d velocity_z_pow_vec = _mm256_mul_pd(velocity_z_vec, velocity_z_vec);        // velocity_field_z[x, y, z] ** 2

        __m256d velocity_x_y_pow_sum_vec = _mm256_add_pd(velocity_x_pow_vec, velocity_y_pow_vec);                       // velocity_field_x[x, y, z] ** 2 + velocity_field_y[x, y, z] ** 2
        __m256d velocity_x_y_z_pow_sum_vec = _mm256_fmadd_pd(velocity_z_vec, velocity_z_vec, velocity_x_y_pow_sum_vec); // velocity_field_x[x, y, z] ** 2 + velocity_field_y[x, y, z] ** 2 + velocity_field_z[x, y, z] ** 2

        __m256d norm_square_vec = _mm256_mul_pd(velocity_x_y_z_pow_sum_vec, const_cs_2_2_inv_vec); // norm_square
        __m256d one_min_norm_square_vec = _mm256_sub_pd(const_1_vec, norm_square_vec);             // 1.0 - norm_square

        __m256d feq_prod_0_vec = _mm256_mul_pd(const_weights_29_vec, density_field_vec); // weights_29 * density_field[x, y, z]
        __m256d feq_0_vec = _mm256_mul_pd(feq_prod_0_vec, one_min_norm_square_vec);      // feq

        __m256d particle_dist_0_vec = _mm256_load_pd(&particle_distributions[scalar_index(x, y, z, 0)]);    // particle_distributions[x, y, z, 0]
        __m256d particle_dist_0_prod_0 = _mm256_mul_pd(const_omtauinv_vec, particle_dist_0_vec);            // omtauinv * particle_distributions[x, y, z, 0]
        __m256d particle_dist_0_res = _mm256_fmadd_pd(const_tauinv_vec, feq_0_vec, particle_dist_0_prod_0); // omtauinv * particle_distributions[x, y, z, 0] + tauinv * feq
        _mm256_store_pd(&particle_distributions[scalar_index(x, y, z, 0)], particle_dist_0_res);

        // ------------- l = 1 -------------

        __m256d vel_x_pow_div_c_s_4_vec = _mm256_mul_pd(velocity_x_pow_vec, const_cs_4_inv_vec);          // velocity_field_x[x, y, z] ** 2 / c_s_4
        __m256d vel_x_sum = _mm256_fmadd_pd(velocity_x_vec, const_cs_2_inv_vec, vel_x_pow_div_c_s_4_vec); // velocity_field_x[x, y, z] ** 2 / c_s_4 + velocity_field_x[x, y, z] / c_s_2

        __m256d particle_dist_1_vec = _mm256_load_pd(&particle_distributions[scalar_index(x, y, z, 1)]); // particle_distributions[x, y, z, 1]
        __m256d particle_dist_1_prod_0 = _mm256_mul_pd(const_omtauinv_vec, particle_dist_1_vec);         // omtauinv * particle_distributions[x, y, z, 1]

        __m256d feq_1_sum = _mm256_add_pd(vel_x_sum, one_min_norm_square_vec); // 1.0 + velocity_field_x[scalar_index(x, y, z)] / c_s_2 + velocity_field_x[scalar_index(x, y, z)] * velocity_field_x[scalar_index(x, y, z)] / c_s_4 - norm_square
        __m256d feq_1_prod = _mm256_mul_pd(feq_1_sum, weights_19_mul_density);
        __m256d feq_1_res = _mm256_fmadd_pd(const_tauinv_vec, feq_1_prod, particle_dist_1_prod_0);
        _mm256_store_pd(&previous_particle_distributions[scalar_index(x, y, z, 1)], feq_1_res);

        // ------------- l = 2 -------------

        __m256d particle_dist_2_vec = _mm256_load_pd(&particle_distributions[scalar_index(x, y, z, 2)]); // particle_distributions[x, y, z, 2]
        __m256d particle_dist_2_prod_0 = _mm256_mul_pd(const_omtauinv_vec, particle_dist_2_vec);         // omtauinv * particle_distributions[x, y, z, 2]

        __m256d vel_x_diff = _mm256_fmadd_pd(velocity_x_vec, const_min_cs_2_inv_vec, vel_x_pow_div_c_s_4_vec); // velocity_field_x[x, y, z] ** 2 / c_s_4 - velocity_field_x[x, y, z] / c_s_2
        __m256d feq_2_sum = _mm256_add_pd(vel_x_diff, one_min_norm_square_vec); // 1.0 + velocity_field_x[scalar_index(x, y, z)] / c_s_2 + velocity_field_x[scalar_index(x, y, z)] * velocity_field_x[scalar_index(x, y, z)] / c_s_4 - norm_square
        __m256d feq_2_prod = _mm256_mul_pd(feq_2_sum, weights_19_mul_density);
        __m256d feq_2_res = _mm256_fmadd_pd(const_tauinv_vec, feq_2_prod, particle_dist_2_prod_0);
        _mm256_store_pd(&previous_particle_distributions[scalar_index(x, y, z, 2)], feq_2_res);

        // ------------- l = 3 -------------

        __m256d vel_y_pow_div_c_s_4_vec = _mm256_mul_pd(velocity_y_pow_vec, const_cs_4_inv_vec);          // velocity_field_y[x, y, z] ** 2 / c_s_4
        __m256d vel_y_sum = _mm256_fmadd_pd(velocity_y_vec, const_cs_2_inv_vec, vel_y_pow_div_c_s_4_vec); // velocity_field_xyx, y, z] ** 2 / c_s_4 + velocity_field_y[x, y, z] / c_s_2

        __m256d particle_dist_3_vec = _mm256_load_pd(&particle_distributions[scalar_index(x, y, z, 3)]); // particle_distributions[x, y, z, 3]
        __m256d particle_dist_3_prod_0 = _mm256_mul_pd(const_omtauinv_vec, particle_dist_3_vec);         // omtauinv * particle_distributions[x, y, z, 3]

        __m256d feq_3_sum = _mm256_add_pd(vel_y_sum, one_min_norm_square_vec); // 1.0 + velocity_field_y[scalar_index(x, y, z)] / c_s_2 + velocity_field_y[scalar_index(x, y, z)] * velocity_field_y[scalar_index(x, y, z)] / c_s_4 - norm_square
        __m256d feq_3_prod = _mm256_mul_pd(feq_3_sum, weights_19_mul_density);
        __m256d feq_3_res = _mm256_fmadd_pd(const_tauinv_vec, feq_3_prod, particle_dist_3_prod_0);
        _mm256_store_pd(&previous_particle_distributions[scalar_index(x, y, z, 3)], feq_3_res);

        // ------------- l = 4 -------------

        __m256d particle_dist_4_vec = _mm256_load_pd(&particle_distributions[scalar_index(x, y, z, 4)]); // particle_distributions[x, y, z, 4]
        __m256d particle_dist_4_prod_0 = _mm256_mul_pd(const_omtauinv_vec, particle_dist_4_vec);         // omtauinv * particle_distributions[x, y, z, 4]

        __m256d vel_y_diff = _mm256_fmadd_pd(velocity_y_vec, const_min_cs_2_inv_vec, vel_y_pow_div_c_s_4_vec); // velocity_field_y[x, y, z] ** 2 / c_s_4 - velocity_field_y[x, y, z] / c_s_2
        __m256d feq_4_sum = _mm256_add_pd(vel_y_diff, one_min_norm_square_vec); // 1.0 + velocity_field_x[scalar_indey(x, y, z)] / c_s_2 + velocity_field_y[scalar_index(x, y, z)] * velocity_field_y[scalar_index(x, y, z)] / c_s_4 - norm_square
        __m256d feq_4_prod = _mm256_mul_pd(feq_4_sum, weights_19_mul_density);
        __m256d feq_4_res = _mm256_fmadd_pd(const_tauinv_vec, feq_4_prod, particle_dist_4_prod_0);
        _mm256_store_pd(&previous_particle_distributions[scalar_index(x, y, z, 4)], feq_4_res);

        // ------------- l = 5 -------------

        __m256d vel_z_pow_div_c_s_4_vec = _mm256_mul_pd(velocity_z_pow_vec, const_cs_4_inv_vec);          // velocity_field_z[x, y, z] ** 2 / c_s_4
        __m256d vel_z_sum = _mm256_fmadd_pd(velocity_z_vec, const_cs_2_inv_vec, vel_z_pow_div_c_s_4_vec); // velocity_field_xzx, y, z] ** 2 / c_s_4 + velocity_field_z[x, y, z] / c_s_2

        __m256d particle_dist_5_vec = _mm256_load_pd(&particle_distributions[scalar_index(x, y, z, 5)]); // particle_distributions[x, y, z, 4]
        __m256d particle_dist_5_prod_0 = _mm256_mul_pd(const_omtauinv_vec, particle_dist_5_vec);         // omtauinv * particle_distributions[x, y, z, 4]

        __m256d feq_5_sum = _mm256_add_pd(vel_z_sum, one_min_norm_square_vec); // 1.0 + velocity_field_z[scalar_index(x, y, z)] / c_s_2 + velocity_field_z[scalar_index(x, y, z)] * velocity_field_y[scalar_index(x, y, z)] / c_s_4 - norm_square
        __m256d feq_5_prod = _mm256_mul_pd(feq_5_sum, weights_19_mul_density);
        __m256d feq_5_res = _mm256_fmadd_pd(const_tauinv_vec, feq_5_prod, particle_dist_5_prod_0);
        _mm256_store_pd(&previous_particle_distributions[scalar_index(x, y, z, 5)], feq_5_res);

        // ------------- l = 6 -------------

        __m256d particle_dist_6_vec = _mm256_load_pd(&particle_distributions[scalar_index(x, y, z, 6)]); // particle_distributions[x, y, z, 6]
        __m256d particle_dist_6_prod_0 = _mm256_mul_pd(const_omtauinv_vec, particle_dist_6_vec);         // omtauinv * particle_distributions[x, y, z, 6]

        __m256d vel_z_diff = _mm256_fmadd_pd(velocity_z_vec, const_min_cs_2_inv_vec, vel_z_pow_div_c_s_4_vec); // velocity_field_y[x, y, z] ** 2 / c_s_6 - velocity_field_y[x, y, z] / c_s_2
        __m256d feq_6_sum = _mm256_add_pd(vel_z_diff, one_min_norm_square_vec); // 1.0 + velocity_field_x[scalar_indey(x, y, z)] / c_s_2 + velocity_field_y[scalar_index(x, y, z)] * velocity_field_y[scalar_index(x, y, z)] / c_s_6 - norm_square
        __m256d feq_6_prod = _mm256_mul_pd(feq_6_sum, weights_19_mul_density);
        __m256d feq_6_res = _mm256_fmadd_pd(const_tauinv_vec, feq_6_prod, particle_dist_6_prod_0);
        _mm256_store_pd(&previous_particle_distributions[scalar_index(x, y, z, 6)], feq_6_res);

        // ------------- l = 7 -------------

        __m256d dot_product_7_sum_0 = _mm256_add_pd(velocity_x_vec, velocity_y_vec);
        __m256d dot_product_7 = _mm256_add_pd(dot_product_7_sum_0, velocity_z_vec);

        __m256d feq_7_prod_0 = _mm256_mul_pd(dot_product_7, dot_product_7);
        __m256d feq_7_prod_1 = _mm256_mul_pd(dot_product_7, const_cs_2_inv_vec);
        __m256d feq_7_sum_0 = _mm256_fmadd_pd(feq_7_prod_0, const_cs_4_inv_vec, one_min_norm_square_vec);
        __m256d feq_7_sum_1 = _mm256_add_pd(feq_7_prod_1, feq_7_sum_0);
        __m256d feq_7 = _mm256_mul_pd(weights_172_mul_density, feq_7_sum_1);

        __m256d particle_dist_7_vec = _mm256_load_pd(&particle_distributions[scalar_index(x, y, z, 7)]); // particle_distributions[x, y, z, 7]
        __m256d particle_dist_7_prod_0 = _mm256_mul_pd(const_omtauinv_vec, particle_dist_7_vec);         // omtauinv * particle_distributions[x, y, z, 7]

        __m256d feq_7_res = _mm256_fmadd_pd(feq_7, const_tauinv_vec, particle_dist_7_prod_0);
        _mm256_store_pd(&previous_particle_distributions[scalar_index(x, y, z, 7)], feq_7_res);

        // ------------- l = 8 -------------

        __m256d feq_8_prod_0 = _mm256_mul_pd(dot_product_7, dot_product_7);
        __m256d feq_8_prod_1 = _mm256_mul_pd(dot_product_7, const_min_cs_2_inv_vec);
        __m256d feq_8_sum_0 = _mm256_fmadd_pd(feq_8_prod_0, const_cs_4_inv_vec, one_min_norm_square_vec);
        __m256d feq_8_sum_1 = _mm256_add_pd(feq_8_prod_1, feq_8_sum_0);
        __m256d feq_8 = _mm256_mul_pd(weights_172_mul_density, feq_8_sum_1);

        __m256d particle_dist_8_vec = _mm256_load_pd(&particle_distributions[scalar_index(x, y, z, 8)]); // particle_distributions[x, y, z, 8]
        __m256d particle_dist_8_prod_0 = _mm256_mul_pd(const_omtauinv_vec, particle_dist_8_vec);         // omtauinv * particle_distributions[x, y, z, 8]

        __m256d feq_8_res = _mm256_fmadd_pd(feq_8, const_tauinv_vec, particle_dist_8_prod_0);
        _mm256_store_pd(&previous_particle_distributions[scalar_index(x, y, z, 8)], feq_8_res);

        // ------------- l = 9 -------------

        __m256d dot_product_9_sum_0 = _mm256_add_pd(velocity_x_vec, velocity_y_vec);
        __m256d dot_product_9 = _mm256_sub_pd(dot_product_9_sum_0, velocity_z_vec);

        __m256d feq_9_prod_0 = _mm256_mul_pd(dot_product_9, dot_product_9);
        __m256d feq_9_prod_1 = _mm256_mul_pd(dot_product_9, const_cs_2_inv_vec);
        __m256d feq_9_sum_0 = _mm256_fmadd_pd(feq_9_prod_0, const_cs_4_inv_vec, one_min_norm_square_vec);
        __m256d feq_9_sum_1 = _mm256_add_pd(feq_9_prod_1, feq_9_sum_0);
        __m256d feq_9 = _mm256_mul_pd(weights_172_mul_density, feq_9_sum_1);

        __m256d particle_dist_9_vec = _mm256_load_pd(&particle_distributions[scalar_index(x, y, z, 9)]); // particle_distributions[x, y, z, 9]
        __m256d particle_dist_9_prod_0 = _mm256_mul_pd(const_omtauinv_vec, particle_dist_9_vec);         // omtauinv * particle_distributions[x, y, z, 9]

        __m256d feq_9_res = _mm256_fmadd_pd(feq_9, const_tauinv_vec, particle_dist_9_prod_0);
        _mm256_store_pd(&previous_particle_distributions[scalar_index(x, y, z, 9)], feq_9_res);

        // ------------- l = 10 -------------

        __m256d feq_10_prod_0 = _mm256_mul_pd(dot_product_9, dot_product_9);
        __m256d feq_10_prod_1 = _mm256_mul_pd(dot_product_9, const_min_cs_2_inv_vec);
        __m256d feq_10_sum_0 = _mm256_fmadd_pd(feq_10_prod_0, const_cs_4_inv_vec, one_min_norm_square_vec);
        __m256d feq_10_sum_1 = _mm256_add_pd(feq_10_prod_1, feq_10_sum_0);
        __m256d feq_10 = _mm256_mul_pd(weights_172_mul_density, feq_10_sum_1);

        __m256d particle_dist_10_vec = _mm256_load_pd(&particle_distributions[scalar_index(x, y, z, 10)]); // particle_distributions[x, y, z, 10]
        __m256d particle_dist_10_prod_0 = _mm256_mul_pd(const_omtauinv_vec, particle_dist_10_vec);         // omtauinv * particle_distributions[x, y, z, 10]

        __m256d feq_10_res = _mm256_fmadd_pd(feq_10, const_tauinv_vec, particle_dist_10_prod_0);
        _mm256_store_pd(&previous_particle_distributions[scalar_index(x, y, z, 10)], feq_10_res);

        // ------------- l = 11 -------------

        __m256d dot_product_11_sum_0 = _mm256_sub_pd(velocity_x_vec, velocity_y_vec);
        __m256d dot_product_11 = _mm256_add_pd(dot_product_11_sum_0, velocity_z_vec);

        __m256d feq_11_prod_0 = _mm256_mul_pd(dot_product_11, dot_product_11);
        __m256d feq_11_prod_1 = _mm256_mul_pd(dot_product_11, const_cs_2_inv_vec);
        __m256d feq_11_sum_0 = _mm256_fmadd_pd(feq_11_prod_0, const_cs_4_inv_vec, one_min_norm_square_vec);
        __m256d feq_11_sum_1 = _mm256_add_pd(feq_11_prod_1, feq_11_sum_0);
        __m256d feq_11 = _mm256_mul_pd(weights_172_mul_density, feq_11_sum_1);

        __m256d particle_dist_11_vec = _mm256_load_pd(&particle_distributions[scalar_index(x, y, z, 11)]); // particle_distributions[x, y, z, 11]
        __m256d particle_dist_11_prod_0 = _mm256_mul_pd(const_omtauinv_vec, particle_dist_11_vec);         // omtauinv * particle_distributions[x, y, z, 11]

        __m256d feq_11_res = _mm256_fmadd_pd(feq_11, const_tauinv_vec, particle_dist_11_prod_0);
        _mm256_store_pd(&previous_particle_distributions[scalar_index(x, y, z, 11)], feq_11_res);

        // ------------- l = 12 -------------

        __m256d feq_12_prod_0 = _mm256_mul_pd(dot_product_11, dot_product_11);
        __m256d feq_12_prod_1 = _mm256_mul_pd(dot_product_11, const_min_cs_2_inv_vec);
        __m256d feq_12_sum_0 = _mm256_fmadd_pd(feq_12_prod_0, const_cs_4_inv_vec, one_min_norm_square_vec);
        __m256d feq_12_sum_1 = _mm256_add_pd(feq_12_prod_1, feq_12_sum_0);
        __m256d feq_12 = _mm256_mul_pd(weights_172_mul_density, feq_12_sum_1);

        __m256d particle_dist_12_vec = _mm256_load_pd(&particle_distributions[scalar_index(x, y, z, 12)]); // particle_distributions[x, y, z, 12]
        __m256d particle_dist_12_prod_0 = _mm256_mul_pd(const_omtauinv_vec, particle_dist_12_vec);         // omtauinv * particle_distributions[x, y, z, 12]

        __m256d feq_12_res = _mm256_fmadd_pd(feq_12, const_tauinv_vec, particle_dist_12_prod_0);
        _mm256_store_pd(&previous_particle_distributions[scalar_index(x, y, z, 12)], feq_12_res);

        // ------------- l = 13 -------------

        __m256d dot_product_13_sum_0 = _mm256_add_pd(velocity_y_vec, velocity_z_vec);
        __m256d dot_product_13 = _mm256_sub_pd(dot_product_13_sum_0, velocity_x_vec);

        __m256d feq_13_prod_0 = _mm256_mul_pd(dot_product_13, dot_product_13);
        __m256d feq_13_prod_1 = _mm256_mul_pd(dot_product_13, const_cs_2_inv_vec);
        __m256d feq_13_sum_0 = _mm256_fmadd_pd(feq_13_prod_0, const_cs_4_inv_vec, one_min_norm_square_vec);
        __m256d feq_13_sum_1 = _mm256_add_pd(feq_13_prod_1, feq_13_sum_0);
        __m256d feq_13 = _mm256_mul_pd(weights_172_mul_density, feq_13_sum_1);

        __m256d particle_dist_13_vec = _mm256_load_pd(&particle_distributions[scalar_index(x, y, z, 13)]); // particle_distributions[x, y, z, 13]
        __m256d particle_dist_13_prod_0 = _mm256_mul_pd(const_omtauinv_vec, particle_dist_13_vec);         // omtauinv * particle_distributions[x, y, z, 13]

        __m256d feq_13_res = _mm256_fmadd_pd(feq_13, const_tauinv_vec, particle_dist_13_prod_0);
        _mm256_store_pd(&previous_particle_distributions[scalar_index(x, y, z, 13)], feq_13_res);

        // ------------- l = 14 -------------

        __m256d feq_14_prod_0 = _mm256_mul_pd(dot_product_13, dot_product_13);
        __m256d feq_14_prod_1 = _mm256_mul_pd(dot_product_13, const_min_cs_2_inv_vec);
        __m256d feq_14_sum_0 = _mm256_fmadd_pd(feq_14_prod_0, const_cs_4_inv_vec, one_min_norm_square_vec);
        __m256d feq_14_sum_1 = _mm256_add_pd(feq_14_prod_1, feq_14_sum_0);
        __m256d feq_14 = _mm256_mul_pd(weights_172_mul_density, feq_14_sum_1);

        __m256d particle_dist_14_vec = _mm256_load_pd(&particle_distributions[scalar_index(x, y, z, 14)]); // particle_distributions[x, y, z, 14]
        __m256d particle_dist_14_prod_0 = _mm256_mul_pd(const_omtauinv_vec, particle_dist_14_vec);         // omtauinv * particle_distributions[x, y, z, 14]

        __m256d feq_14_res = _mm256_fmadd_pd(feq_14, const_tauinv_vec, particle_dist_14_prod_0);
        _mm256_store_pd(&previous_particle_distributions[scalar_index(x, y, z, 14)], feq_14_res);
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
        double percentage = (double)(i + 1) / (double)(runs) * 100.0;
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
#include <iostream>
// Now Linux only.
#include "profile.h"
#include "utils.h"
#include <fstream>
#include <iostream>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

// the variables are in utils.h
using namespace std;

// for now one can define OUTPUT, TIMING
// #define OUTPUT

// #define TIMING
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

// Direction Size
#define DS 15

double *density_field;
double *velocity_field;
double *particle_distributions;
double *prev_3;
double *prev_7;
double *prev_9;
double *prev_12;
double *prev_13;

inline int scalar_index(int x, int y, int z) { return (z * NX * NY) + (y * NX) + x; }

inline int scalar_index(int x, int y, int z, int w) { return (x + y * NX + z * NX * NY + w * NX * NY * NZ); }

void initialise() {
  int box_flatten_length = NX * NY * NZ;
  int distributions_flatten_length = box_flatten_length * DS;
  density_field = (double *)malloc(box_flatten_length * sizeof(double));
  velocity_field = (double *)malloc(box_flatten_length * 3 * sizeof(double));
  particle_distributions = (double *)malloc(distributions_flatten_length * sizeof(double));

  prev_3 = (double *)malloc(NX * sizeof(double));
  prev_7 = (double *)malloc(NX * sizeof(double));
  prev_9 = (double *)malloc(NX * sizeof(double));
  prev_12 = (double *)malloc(NX * sizeof(double));
  prev_13 = (double *)malloc(NX * sizeof(double));

  for (int i = 0; i < NX * NY * NZ; i++) {
    density_field[i] = 1.0;
  }

  for (int z = 0; z < NZ; z++) {
    for (int y = 0; y < NY; y++) {
      for (int x = 0; x < NX; x++) {

        for (int i = 0; i < DS; i++) {
          //   previous_particle_distributions[scalar_index(x, y, z, i)] = weights[i];
          if (i == 0)
            particle_distributions[scalar_index(x, y, z, i)] = weights_29;
          else if (i < 7)
            particle_distributions[scalar_index(x, y, z, i)] = weights_19;
          else
            particle_distributions[scalar_index(x, y, z, i)] = weights_172;
        }
      }
    }
  }

  compute_density_momentum_profiler = init_profiler(NZ * NY * NX * (14 * 4 + 2), NX * NY * NZ * 8 * 19);
  collision_profiler = init_profiler(NZ * NY * NX * 130, NX * NY * NZ * 8 * (3 + 3 + 15 + 15) + 15 * 8);
  stream_profiler = init_profiler(NZ * 8, 15 * 2 * NZ * NY * NX + NZ * 10 + NZ * 18);
}

#define c_s_4 (2 * c_s * c_s * c_s * c_s)
#define c_s_2 (c_s * c_s)
#define c_s_2_inv (1 / c_s_2)
#define c_s_4_inv (1 / c_s_4)
#define c_s_2_inv_2 (1 / (2 * c_s_2))
#define c_s_2_inv_02 (c_s_2_inv * 0.2)

int box_length = NX * NY * NZ;
void compute_density_momentum_moment() {
  // negation doesnt count as a flop, right?
  // flops = NZ*NY*NX*(14*4+2)
  // bytes = NX*NY*NZ*8*19
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

        double dens_inv = 1 / new_density;
        density_field[scalar_index(x, y, z)] = new_density;
        velocity_field[(scalar_index(x, y, z)) * 3] = x_sum * dens_inv;
        velocity_field[(scalar_index(x, y, z)) * 3 + 1] = y_sum * dens_inv;
        velocity_field[(scalar_index(x, y, z)) * 3 + 2] = z_sum * dens_inv;
      }
    }
  }
}

void stream() {}

#define weights_29_tau (weights_29 * tauinv)
#define weights_19_tau (weights_19 * tauinv)
#define weights_172_tau (weights_172 * tauinv)

void collision() {
  for (int z = 0; z < NZ; z++) {
    int y = 0;
    for (int x = 0; x < NX; x++) {
      // 13
      int index = scalar_index(x, y, z);
      double norm_square = velocity_field[index * 3] * velocity_field[index * 3] + velocity_field[index * 3 + 1] * velocity_field[index * 3 + 1] + velocity_field[index * 3 + 2] * velocity_field[index * 3 + 2];
      double norm_square_cs2 = 1.0 - norm_square * c_s_2_inv_2;
      double density_field_29 = density_field[index] * weights_29_tau;
      double density_field_19 = density_field[index] * weights_19_tau;
      double density_field_172 = density_field[index] * weights_172_tau;
      prev_3[x] = weights_19;
      prev_7[x] = weights_172;
      prev_12[x] = weights_172;
      prev_13[x] = weights_172;
      prev_9[x] = weights_172;

      double feq_0 = density_field_29 * norm_square_cs2;
      particle_distributions[scalar_index(x, y, z, 0)] = omtauinv * particle_distributions[scalar_index(x, y, z, 0)] + feq_0;

      double feq_1 = density_field_19 * (velocity_field[index * 3] * (c_s_2_inv + velocity_field[index * 3] * c_s_4_inv) + norm_square_cs2);
      particle_distributions[scalar_index(x, y, z, 1)] = omtauinv * particle_distributions[scalar_index(x, y, z, 1)] + feq_1;

      double dot_product_2 = -velocity_field[index * 3];
      double feq_2 = density_field_19 * (dot_product_2 * (c_s_2_inv + dot_product_2 * c_s_4_inv) + norm_square_cs2);
      particle_distributions[scalar_index(x, y, z, 2)] = omtauinv * particle_distributions[scalar_index(x, y, z, 2)] + feq_2;
      // 7
      double dot_product_5 = velocity_field[index * 3 + 2];
      double feq_5 = density_field_19 * (dot_product_5 * (c_s_2_inv + dot_product_5 * c_s_4_inv) + norm_square_cs2);
      particle_distributions[scalar_index(x, y, z, 5)] = omtauinv * particle_distributions[scalar_index(x, y, z, 5)] + feq_5;

      // 7
      double dot_product_6 = velocity_field[index * 3 + 2];
      double feq_6 = density_field_19 * (dot_product_6 * (-c_s_2_inv + dot_product_6 * c_s_4_inv) + norm_square_cs2);
      particle_distributions[scalar_index(x, y, z, 6)] = omtauinv * particle_distributions[scalar_index(x, y, z, 6)] + feq_6;
    }

    for (int y = 1; y < NY - 1; y++) {
      for (int x = 0; x < NX; x++) {
        // 13
        int index = scalar_index(x, y, z);
        double norm_square = velocity_field[index * 3] * velocity_field[index * 3] + velocity_field[index * 3 + 1] * velocity_field[index * 3 + 1] + velocity_field[index * 3 + 2] * velocity_field[index * 3 + 2];
        double norm_square_cs2 = 1.0 - norm_square * c_s_2_inv_2;
        double density_field_29 = density_field[index] * weights_29_tau;
        double density_field_19 = density_field[index] * weights_19_tau;
        double density_field_172 = density_field[index] * weights_172_tau;

        // 3
        double feq_0 = density_field_29 * norm_square_cs2;
        particle_distributions[scalar_index(x, y, z, 0)] = omtauinv * particle_distributions[scalar_index(x, y, z, 0)] + feq_0;
        // 7
        double feq_1 = density_field_19 * (velocity_field[index * 3] * (c_s_2_inv + velocity_field[index * 3] * c_s_4_inv) + norm_square_cs2);
        particle_distributions[scalar_index(x, y, z, 1)] = omtauinv * particle_distributions[scalar_index(x, y, z, 1)] + feq_1;

        // 7
        double feq_2 = density_field_19 * (velocity_field[index * 3] * (-c_s_2_inv + velocity_field[index * 3] * c_s_4_inv) + norm_square_cs2);
        particle_distributions[scalar_index(x, y, z, 2)] = omtauinv * particle_distributions[scalar_index(x, y, z, 2)] + feq_2;

        // 7
        double feq_3 = density_field_19 * (velocity_field[index * 3 + 1] * (c_s_2_inv + velocity_field[index * 3 + 1] * c_s_4_inv) + norm_square_cs2);
        double curr_3 = omtauinv * particle_distributions[scalar_index(x, y, z, 3)] + feq_3;
        particle_distributions[scalar_index(x, y, z, 3)] = prev_3[x];
        prev_3[x] = curr_3;

        // 7
        double dot_product_4 = -velocity_field[index * 3 + 1];
        double feq_4 = density_field_19 * (dot_product_4 * (c_s_2_inv + dot_product_4 * c_s_4_inv) + norm_square_cs2);
        particle_distributions[scalar_index(x, y - 1, z, 4)] = omtauinv * particle_distributions[scalar_index(x, y, z, 4)] + feq_4;

        // 7
        double dot_product_5 = velocity_field[index * 3 + 2];
        double feq_5 = density_field_19 * (dot_product_5 * (c_s_2_inv + dot_product_5 * c_s_4_inv) + norm_square_cs2);
        particle_distributions[scalar_index(x, y, z, 5)] = omtauinv * particle_distributions[scalar_index(x, y, z, 5)] + feq_5;

        // 7
        double dot_product_6 = -velocity_field[index * 3 + 2];
        double feq_6 = density_field_19 * (dot_product_6 * (c_s_2_inv + dot_product_6 * c_s_4_inv) + norm_square_cs2);
        particle_distributions[scalar_index(x, y, z, 6)] = omtauinv * particle_distributions[scalar_index(x, y, z, 6)] + feq_6;

        // 9
        double dot_product_7 = velocity_field[index * 3] + velocity_field[index * 3 + 1] + velocity_field[index * 3 + 2];
        double feq_7 = density_field_172 * (dot_product_7 * (c_s_2_inv + dot_product_7 * c_s_4_inv) + norm_square_cs2);

        double curr_7 = omtauinv * particle_distributions[scalar_index(x, y, z, 7)] + feq_7;
        particle_distributions[scalar_index(x, y, z, 7)] = prev_7[x];
        prev_7[x] = curr_7;

        // 9
        double dot_product_8 = -velocity_field[index * 3] - velocity_field[index * 3 + 1] - velocity_field[index * 3 + 2];
        double feq_8 = density_field_172 * (dot_product_8 * (c_s_2_inv + dot_product_8 * c_s_4_inv) + norm_square_cs2);
        // TODO: this type of thing can def be rewrtitten somehow for sure
        particle_distributions[scalar_index(x, y - 1, z, 8)] = omtauinv * particle_distributions[scalar_index(x, y, z, 8)] + feq_8;

        // 9
        double dot_product_9 = velocity_field[index * 3] + velocity_field[index * 3 + 1] - velocity_field[index * 3 + 2];
        double feq_9 = density_field_172 * (dot_product_9 * (c_s_2_inv + dot_product_9 * c_s_4_inv) + norm_square_cs2);

        double curr_9 = omtauinv * particle_distributions[scalar_index(x, y, z, 9)] + feq_9;
        particle_distributions[scalar_index(x, y, z, 9)] = prev_9[x];
        prev_9[x] = curr_9;

        // 9
        double dot_product_10 = -velocity_field[index * 3] - velocity_field[index * 3 + 1] + velocity_field[index * 3 + 2];
        double feq_10 = density_field_172 * (dot_product_10 * (c_s_2_inv + dot_product_10 * c_s_4_inv) + norm_square_cs2);
        particle_distributions[scalar_index(x, y - 1, z, 10)] = omtauinv * particle_distributions[scalar_index(x, y, z, 10)] + feq_10;

        // 9
        double dot_product_11 = velocity_field[index * 3] - velocity_field[index * 3 + 1] + velocity_field[index * 3 + 2];
        double feq_11 = density_field_172 * (dot_product_11 * (c_s_2_inv + dot_product_11 * c_s_4_inv) + norm_square_cs2);
        particle_distributions[scalar_index(x, y - 1, z, 11)] = omtauinv * particle_distributions[scalar_index(x, y, z, 11)] + feq_11;

        // 9
        double dot_product_12 = -velocity_field[index * 3] + velocity_field[index * 3 + 1] - velocity_field[index * 3 + 2];
        double feq_12 = density_field_172 * (dot_product_12 * (c_s_2_inv + dot_product_12 * c_s_4_inv) + norm_square_cs2);
        double curr_12 = omtauinv * particle_distributions[scalar_index(x, y, z, 12)] + feq_12;
        particle_distributions[scalar_index(x, y, z, 12)] = prev_12[x];
        prev_12[x] = curr_12;

        // 9
        double dot_product_13 = -velocity_field[index * 3] + velocity_field[index * 3 + 1] + velocity_field[index * 3 + 2];
        double feq_13 = density_field_172 * (dot_product_13 * (c_s_2_inv + dot_product_13 * c_s_4_inv) + norm_square_cs2);
        double curr_13 = omtauinv * particle_distributions[scalar_index(x, y, z, 13)] + feq_13;
        particle_distributions[scalar_index(x, y, z, 13)] = prev_13[x];
        prev_13[x] = curr_13;

        // 9
        double dot_product_14 = velocity_field[index * 3] - velocity_field[index * 3 + 1] - velocity_field[index * 3 + 2];
        double feq_14 = density_field_172 * (dot_product_14 * (c_s_2_inv + dot_product_14 * c_s_4_inv) + norm_square_cs2);
        particle_distributions[scalar_index(x, y - 1, z, 14)] = omtauinv * particle_distributions[scalar_index(x, y, z, 14)] + feq_14;
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

        double dens_inv = 1 / new_density;
        density_field[scalar_index(x, y, z)] = new_density;
        velocity_field[(scalar_index(x, y, z)) * 3] = x_sum * dens_inv;
        velocity_field[(scalar_index(x, y, z)) * 3 + 1] = y_sum * dens_inv;
        velocity_field[(scalar_index(x, y, z)) * 3 + 2] = z_sum * dens_inv;
      }
    }
    y = NY - 1;
    for (int x = 0; x < NX; x++) {
      // 13
      int index = scalar_index(x, y, z);
      double norm_square = velocity_field[index * 3] * velocity_field[index * 3] + velocity_field[index * 3 + 1] * velocity_field[index * 3 + 1] + velocity_field[index * 3 + 2] * velocity_field[index * 3 + 2];
      double norm_square_cs2 = 1.0 - norm_square * c_s_2_inv_2;
      double density_field_29 = density_field[index] * weights_29_tau;
      double density_field_19 = density_field[index] * weights_19_tau;
      double density_field_172 = density_field[index] * weights_172_tau;

      // 3
      double feq_0 = density_field_29 * norm_square_cs2;
      particle_distributions[scalar_index(x, y, z, 0)] = omtauinv * particle_distributions[scalar_index(x, y, z, 0)] + feq_0;
      // 7
      double dot_product_1 = velocity_field[index * 3];
      double feq_1 = density_field_19 * (dot_product_1 * (c_s_2_inv + dot_product_1 * c_s_4_inv) + norm_square_cs2);
      particle_distributions[scalar_index(x, y, z, 1)] = omtauinv * particle_distributions[scalar_index(x, y, z, 1)] + feq_1;

      // 7
      double dot_product_2 = -velocity_field[index * 3];
      double feq_2 = density_field_19 * (dot_product_2 * (c_s_2_inv + dot_product_2 * c_s_4_inv) + norm_square_cs2);
      particle_distributions[scalar_index(x, y, z, 2)] = omtauinv * particle_distributions[scalar_index(x, y, z, 2)] + feq_2;

      // 7
      double dot_product_3 = velocity_field[index * 3 + 1];
      double feq_3 = density_field_19 * (dot_product_3 * (c_s_2_inv + dot_product_3 * c_s_4_inv) + norm_square_cs2);
      particle_distributions[scalar_index(x, y, z, 3)] = prev_3[x];

      // 7
      double dot_product_4 = -velocity_field[index * 3 + 1];
      double feq_4 = density_field_19 * (dot_product_4 * (c_s_2_inv + dot_product_4 * c_s_4_inv) + norm_square_cs2);
      particle_distributions[scalar_index(x, y - 1, z, 4)] = omtauinv * particle_distributions[scalar_index(x, y, z, 4)] + feq_4;
      particle_distributions[scalar_index(x, y, z, 4)] = omtauinv * particle_distributions[scalar_index(x, y, z, 3)] + feq_3;

      // 7
      double dot_product_5 = velocity_field[index * 3 + 2];
      double feq_5 = density_field_19 * (dot_product_5 * (c_s_2_inv + dot_product_5 * c_s_4_inv) + norm_square_cs2);
      particle_distributions[scalar_index(x, y, z, 5)] = omtauinv * particle_distributions[scalar_index(x, y, z, 5)] + feq_5;

      // 7
      double dot_product_6 = -velocity_field[index * 3 + 2];
      double feq_6 = density_field_19 * (dot_product_6 * (c_s_2_inv + dot_product_6 * c_s_4_inv) + norm_square_cs2);
      particle_distributions[scalar_index(x, y, z, 6)] = omtauinv * particle_distributions[scalar_index(x, y, z, 6)] + feq_6;

      // 9
      double dot_product_7 = velocity_field[index * 3] + velocity_field[index * 3 + 1] + velocity_field[index * 3 + 2];
      double feq_7 = density_field_172 * (dot_product_7 * (c_s_2_inv + dot_product_7 * c_s_4_inv) + norm_square_cs2);

      double curr_7 = omtauinv * particle_distributions[scalar_index(x, y, z, 7)] + feq_7;
      particle_distributions[scalar_index(x, y, z, 7)] = prev_7[x];
      prev_7[x] = curr_7;
      // 9
      double dot_product_8 = -velocity_field[index * 3] - velocity_field[index * 3 + 1] - velocity_field[index * 3 + 2];
      double feq_8 = density_field_172 * (dot_product_8 * (c_s_2_inv + dot_product_8 * c_s_4_inv) + norm_square_cs2);
      particle_distributions[scalar_index(x, y - 1, z, 8)] = omtauinv * particle_distributions[scalar_index(x, y, z, 8)] + feq_8;
      particle_distributions[scalar_index(x, y, z, 8)] = omtauinv * particle_distributions[scalar_index(x, y, z, 7)] + feq_7 - weights_172 * c_s_2_inv_02;

      // 9
      double dot_product_9 = velocity_field[index * 3] + velocity_field[index * 3 + 1] - velocity_field[index * 3 + 2];
      double feq_9 = density_field_172 * (dot_product_9 * (c_s_2_inv + dot_product_9 * c_s_4_inv) + norm_square_cs2);
      particle_distributions[scalar_index(x, y, z, 9)] = prev_9[x];

      // 9
      double dot_product_10 = -velocity_field[index * 3] - velocity_field[index * 3 + 1] + velocity_field[index * 3 + 2];
      double feq_10 = density_field_172 * (dot_product_10 * (c_s_2_inv + dot_product_10 * c_s_4_inv) + norm_square_cs2);
      particle_distributions[scalar_index(x, y - 1, z, 10)] = omtauinv * particle_distributions[scalar_index(x, y, z, 10)] + feq_10;
      particle_distributions[scalar_index(x, y, z, 10)] = omtauinv * particle_distributions[scalar_index(x, y, z, 9)] + feq_9 - weights_172 * c_s_2_inv_02;

      // 9
      double dot_product_11 = velocity_field[index * 3] - velocity_field[index * 3 + 1] + velocity_field[index * 3 + 2];
      double feq_11 = density_field_172 * (dot_product_11 * (c_s_2_inv + dot_product_11 * c_s_4_inv) + norm_square_cs2);
      particle_distributions[scalar_index(x, y - 1, z, 11)] = omtauinv * particle_distributions[scalar_index(x, y, z, 11)] + feq_11;

      // 9
      double dot_product_12 = -velocity_field[index * 3] + velocity_field[index * 3 + 1] - velocity_field[index * 3 + 2];
      double feq_12 = density_field_172 * (dot_product_12 * (c_s_2_inv + dot_product_12 * c_s_4_inv) + norm_square_cs2);
      particle_distributions[scalar_index(x, y, z, 12)] = prev_12[x];
      particle_distributions[scalar_index(x, y, z, 11)] = omtauinv * particle_distributions[scalar_index(x, y, z, 12)] + feq_12 + weights_172 * c_s_2_inv_02;
      // 9
      double dot_product_13 = -velocity_field[index * 3] + velocity_field[index * 3 + 1] + velocity_field[index * 3 + 2];
      double feq_13 = density_field_172 * (dot_product_13 * (c_s_2_inv + dot_product_13 * c_s_4_inv) + norm_square_cs2);
      particle_distributions[scalar_index(x, y, z, 13)] = prev_13[x];

      // 9
      double dot_product_14 = velocity_field[index * 3] - velocity_field[index * 3 + 1] - velocity_field[index * 3 + 2];
      double feq_14 = density_field_172 * (dot_product_14 * (c_s_2_inv + dot_product_14 * c_s_4_inv) + norm_square_cs2);
      particle_distributions[scalar_index(x, y - 1, z, 14)] = omtauinv * particle_distributions[scalar_index(x, y, z, 14)] + feq_14;
      particle_distributions[scalar_index(x, y, z, 14)] = omtauinv * particle_distributions[scalar_index(x, y, z, 13)] + feq_13 + weights_172 * c_s_2_inv_02;
    }
  }
}

void perform_timestep() {
  time_lbm++;

  start_run(compute_density_momentum_profiler);
  compute_density_momentum_moment();
  end_run(compute_density_momentum_profiler);

  start_run(collision_profiler);
  collision();
  end_run(collision_profiler);

  // start_run(stream_profiler);
  // stream();
  // end_run(stream_profiler);
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
         stream_stats.performance, 0, stream_stats.runs, stream_stats.arithmetic_intensity);
#endif
  return 0;
}
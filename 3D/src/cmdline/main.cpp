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
double *previous_particle_distributions;
double *particle_distributions;

inline int scalar_index(int x, int y, int z) { return (z * NX * NY) + (y * NX) + x; }
inline int scalar_index(int x, int y, int z, int w) { return (x + y * NX + z * NX * NY + w * NX * NY * NZ); }
// inline int scalar_index(int x, int y, int z, int w) { return (z * NX * NY * DS) + (y * NX * DS) + (x * DS) + w; }

// Lattice directions using D3DQ15. assumed speed of sound c_s = 1/sqrt(3).
vector_3_int *directions;
double *weights;
int *reverse_indexes;

void initialise() {
  int box_flatten_length = NX * NY * NZ;
  int distributions_flatten_length = box_flatten_length * DS;
  density_field = (double *)malloc(box_flatten_length * sizeof(double));
  velocity_field = (double *)malloc(box_flatten_length * 3 * sizeof(double));

  previous_particle_distributions = (double *)malloc(distributions_flatten_length * sizeof(double));
  particle_distributions = (double *)malloc(distributions_flatten_length * sizeof(double));

  for (int i = 0; i < NX * NY * NZ; i++) {
    density_field[i] = 1;
  }
  for (int z = 0; z < NZ; z++) {
    for (int y = 0; y < NY; y++) {
      for (int x = 0; x < NX; x++) {

        for (int i = 0; i < DS; i++) {
          previous_particle_distributions[scalar_index(x, y, z, i)] = weights[i];
          particle_distributions[scalar_index(x, y, z, i)] = weights[i];
        }
      }
    }
  }

  compute_density_momentum_profiler = init_profiler(NZ * NY * NX * (14 * 4 + 2), NX * NY * NZ * 8 * 19);
  collision_profiler = init_profiler(NZ * NY * NX * 130, NX * NY * NZ * 8 * (3 + 3 + 15 + 15) + 15 * 8);
  stream_profiler = init_profiler(NZ * 8, 15 * 2 * NZ * NY * NX + NZ * 10 + NZ * 18);
}

#define c_s_4 (2*c_s*c_s*c_s*c_s)
#define c_s_2 (c_s*c_s)
#define c_s_2_inv (1/c_s_2)
#define c_s_4_inv (1/c_s_4)
#define c_s_2_inv_2 (1/(2 * c_s_2))
#define c_s_2_inv_02 (c_s_2_inv*0.2)

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

int box_length = NX*NY*NZ;
void compute_density_momentum_moment() {
  // negation doesnt count as a flop, right?
  // flops = NZ*NY*NX*(14*4+2)
  // bytes = NX*NY*NZ*8*19
  int scalar_index_curr_0, scalar_index_curr_1, scalar_index_curr_2, scalar_index_curr_3, scalar_index_curr_4, scalar_index_curr_5, scalar_index_curr_6, scalar_index_curr_7, scalar_index_curr_8, scalar_index_curr_9, scalar_index_curr_10, scalar_index_curr_11, scalar_index_curr_12, scalar_index_curr_13, scalar_index_curr_14;
  for (int z = 0; z < NZ; z++) {
    for (int y = 0; y < NY; y++) {
      for (int x = 0; x < NX; x++) {
        //we touch 17 blocks here, all entries in order and once.. should be fine?
        double new_density, new_density_0 = 0, new_density_1 = 0, new_density_2 = 0, new_density_3 = 0, new_density_4 = 0, new_density_5 = 0, new_density_6 = 0, new_density_7 = 0, new_density_8 = 0, new_density_9 = 0, new_density_10 = 0,
                            new_density_11 = 0, new_density_12 = 0, new_density_13 = 0, new_density_14 = 0;

        scalar_index_curr_0 = scalar_index(x, y, z, 0);
        scalar_index_curr_1 = scalar_index(x, y, z, 1);
        scalar_index_curr_2 = scalar_index(x, y, z, 2);
        scalar_index_curr_3 = scalar_index(x, y, z, 3);
        scalar_index_curr_4 = scalar_index(x, y, z, 4);
        scalar_index_curr_5 = scalar_index(x, y, z, 5);
        scalar_index_curr_6 = scalar_index(x, y, z, 6);
        scalar_index_curr_7 = scalar_index(x, y, z, 7);
        scalar_index_curr_8 = scalar_index(x, y, z, 8);
        scalar_index_curr_9 = scalar_index(x, y, z, 9);
        scalar_index_curr_10 = scalar_index(x, y, z, 10);
        scalar_index_curr_11 = scalar_index(x, y, z, 11);
        scalar_index_curr_12 = scalar_index(x, y, z, 12);
        scalar_index_curr_13 = scalar_index(x, y, z, 13);
        scalar_index_curr_14 = scalar_index(x, y, z, 14);
        new_density_0 = particle_distributions[scalar_index_curr_0];

        new_density_1 = particle_distributions[scalar_index_curr_1];
        double x_1 = particle_distributions[scalar_index_curr_1];

        new_density_2 = particle_distributions[scalar_index_curr_2];
        double x_2 = -particle_distributions[scalar_index_curr_2];

        new_density_3 = particle_distributions[scalar_index_curr_3];
        double y_3 = particle_distributions[scalar_index_curr_3];

        new_density_4 = particle_distributions[scalar_index_curr_4];
        double y_4 = -particle_distributions[scalar_index_curr_4];

        new_density_5 = particle_distributions[scalar_index_curr_5];
        double z_5 = particle_distributions[scalar_index_curr_5];

        new_density_6 = particle_distributions[scalar_index_curr_6];
        double z_6 = -particle_distributions[scalar_index_curr_6];

        new_density_7 = particle_distributions[scalar_index_curr_7];
        double x_7 = particle_distributions[scalar_index_curr_7];
        double y_7 = particle_distributions[scalar_index_curr_7];
        double z_7 = particle_distributions[scalar_index_curr_7];

        new_density_8 = particle_distributions[scalar_index_curr_8];
        double x_8 = -particle_distributions[scalar_index_curr_8];
        double y_8 = -particle_distributions[scalar_index_curr_8];
        double z_8 = -particle_distributions[scalar_index_curr_8];

        new_density_9 = particle_distributions[scalar_index_curr_9];
        double x_9 = particle_distributions[scalar_index_curr_9];
        double y_9 = particle_distributions[scalar_index_curr_9];
        double z_9 = -particle_distributions[scalar_index_curr_9];

        new_density_10 = particle_distributions[scalar_index_curr_10];
        double x_10 = -particle_distributions[scalar_index_curr_10];
        double y_10 = -particle_distributions[scalar_index_curr_10];
        double z_10 = particle_distributions[scalar_index_curr_10];

        new_density_11 = particle_distributions[scalar_index_curr_11];
        double x_11 = particle_distributions[scalar_index_curr_11];
        double y_11 = -particle_distributions[scalar_index_curr_11];
        double z_11 = particle_distributions[scalar_index_curr_11];

        new_density_12 = particle_distributions[scalar_index_curr_12];
        double x_12 = -particle_distributions[scalar_index_curr_12];
        double y_12 = particle_distributions[scalar_index_curr_12];
        double z_12 = -particle_distributions[scalar_index_curr_12];

        new_density_13 = particle_distributions[scalar_index_curr_13];
        double x_13 = -particle_distributions[scalar_index_curr_13];
        double y_13 = particle_distributions[scalar_index_curr_13];
        double z_13 = particle_distributions[scalar_index_curr_13];

        new_density_14 = particle_distributions[scalar_index_curr_14];
        double x_14 = particle_distributions[scalar_index_curr_14];
        double y_14 = -particle_distributions[scalar_index_curr_14];
        double z_14 = -particle_distributions[scalar_index_curr_14];

        new_density =
            new_density_0 + new_density_1 + new_density_2 + new_density_3 + new_density_4 + new_density_5 + new_density_6 + new_density_7 + new_density_8 + new_density_9 + new_density_10 + new_density_11 + new_density_12 + new_density_13 + new_density_14;
        double x_sum = x_1 + x_2 + x_7 + x_8 + x_9 + x_10 + x_11 + x_12 + x_13 + x_14;
        double y_sum = y_3 + y_4 + y_7 + y_8 + y_9 + y_10 + y_11 + y_12 + y_13 + y_14;
        double z_sum = z_5 + z_6 + z_7 + z_8 + z_9 + z_10 + z_11 + z_12 + z_13 + z_14;

        double dens_inv = 1/new_density;
        density_field[scalar_index(x, y, z)] = new_density;
        velocity_field[(scalar_index(x, y, z))*3] = x_sum * dens_inv;
        velocity_field[(scalar_index(x, y, z))*3+1] = y_sum * dens_inv;
        velocity_field[(scalar_index(x, y, z))*3+2] = z_sum * dens_inv;
      }      
    }
  }
}

void stream() {
  // we can do better memory accesses here most probably...
  for (int z = 0; z < NZ; z++) {
    const int zp1 = (z + 1) % NZ;
    const int zm1 = (NZ + z - 1) % NZ;
    for (int y = 0; y < NY; y++) {
      for (int x = 0; x < NX; x++) {
        const int xp1 = (x + 1) % NX;
        // flops = NZ*8
        // bytes = 15*2*NZ*NY*NX+NZ*10+NZ*18
        particle_distributions[scalar_index(x, y, z, 0)] = previous_particle_distributions[scalar_index(x, y, z, 0)];
        particle_distributions[scalar_index(x, y, z, 1)] = previous_particle_distributions[scalar_index(x, y, z, 1)];
        particle_distributions[scalar_index(x, y, z, 2)] = previous_particle_distributions[scalar_index(xp1, y, z, 2)];
        particle_distributions[scalar_index(x, y, z, 3)] = previous_particle_distributions[scalar_index(x, y - 1, z, 3)];
        particle_distributions[scalar_index(x, y, z, 4)] = previous_particle_distributions[scalar_index(x, y + 1, z, 4)];
        particle_distributions[scalar_index(x, y, z, 5)] = previous_particle_distributions[scalar_index(x, y, zm1, 5)];
        particle_distributions[scalar_index(x, y, z, 6)] = previous_particle_distributions[scalar_index(x, y, zp1, 6)];
        particle_distributions[scalar_index(x, y, z, 7)] = previous_particle_distributions[scalar_index(x, y - 1, zm1, 7)];
        particle_distributions[scalar_index(x, y, z, 8)] = previous_particle_distributions[scalar_index(xp1, y + 1, zp1, 8)];
        particle_distributions[scalar_index(x, y, z, 9)] = previous_particle_distributions[scalar_index(x, y - 1, zp1, 9)];
        particle_distributions[scalar_index(x, y, z, 10)] = previous_particle_distributions[scalar_index(xp1, y + 1, zm1, 10)];
        particle_distributions[scalar_index(x, y, z, 11)] = previous_particle_distributions[scalar_index(x, y + 1, zm1, 11)];
        particle_distributions[scalar_index(x, y, z, 12)] = previous_particle_distributions[scalar_index(xp1, y - 1, zp1, 12)];
        particle_distributions[scalar_index(x, y, z, 13)] = previous_particle_distributions[scalar_index(xp1, y - 1, zm1, 13)];
        particle_distributions[scalar_index(x, y, z, 14)] = previous_particle_distributions[scalar_index(x, y + 1, zp1, 14)];

        if (y == 0) {
          particle_distributions[scalar_index(x, y, z, 3)] = previous_particle_distributions[scalar_index(x, y, z, 4)];
          particle_distributions[scalar_index(x, y, z, 7)] = previous_particle_distributions[scalar_index(x, y, z, 8)];
          particle_distributions[scalar_index(x, y, z, 9)] = previous_particle_distributions[scalar_index(x, y, z, 10)];
          particle_distributions[scalar_index(x, y, z, 12)] = previous_particle_distributions[scalar_index(x, y, z, 11)];
          particle_distributions[scalar_index(x, y, z, 13)] = previous_particle_distributions[scalar_index(x, y, z, 14)];
        }

        if (y == NY - 1) {
          particle_distributions[scalar_index(x, y, z, 4)] = previous_particle_distributions[scalar_index(x, y, z, 3)];
          particle_distributions[scalar_index(x, y, z, 8)] = previous_particle_distributions[scalar_index(x, y, z, 7)] - weights_172 * c_s_2_inv_02;
          particle_distributions[scalar_index(x, y, z, 10)] = previous_particle_distributions[scalar_index(x, y, z, 9)] - weights_172 * c_s_2_inv_02;
          particle_distributions[scalar_index(x, y, z, 11)] = previous_particle_distributions[scalar_index(x, y, z, 12)] + weights_172 * c_s_2_inv_02;
          particle_distributions[scalar_index(x, y, z, 14)] = previous_particle_distributions[scalar_index(x, y, z, 13)] + weights_172 * c_s_2_inv_02;
        }
      }
    }
  }
}

void collision() { 
  // flops = NZ*NY*NX*130
  // bytes = NX*NY*NZ*8*(3+ 3+15+15) + 15*8
  for (int z = 0; z < NZ; z++) {
    for (int y = 0; y < NY; y++) {
      for (int x = 0; x < NX; x++) {
        // 13
        int index = scalar_index(x, y, z);
        double norm_square = velocity_field[index*3] * velocity_field[index*3] + velocity_field[index*3+1] * velocity_field[index*3+1] + velocity_field[index*3+2] * velocity_field[index*3+2];
        double norm_square_cs2 = 1.0 - norm_square * c_s_2_inv_2;
        double density_field_29 = density_field[index]*weights_29*tauinv;
        double density_field_19 = density_field[index]*weights_19*tauinv;
        double density_field_172 = density_field[index]*weights_172*tauinv;

            // 3
            double feq_0 = density_field_29 * norm_square_cs2;
            previous_particle_distributions[scalar_index(x, y, z, 0)] = omtauinv * particle_distributions[scalar_index(x, y, z, 0)] + feq_0;

            // 7
        double dot_product_1 = velocity_field[index*3];
        double feq_1 = density_field_19 * (dot_product_1 * (c_s_2_inv + dot_product_1 * c_s_4_inv) + norm_square_cs2);
        previous_particle_distributions[scalar_index(x, y, z, 1)] = omtauinv * particle_distributions[scalar_index(x, y, z, 1)] + feq_1;

            // 7
        double dot_product_2 = -velocity_field[index*3];
        double feq_2 = density_field_19 * (dot_product_2 * (c_s_2_inv + dot_product_2 * c_s_4_inv) + norm_square_cs2);
        previous_particle_distributions[scalar_index(x, y, z, 2)] = omtauinv * particle_distributions[scalar_index(x, y, z, 2)] + feq_2;

        // 7
        double dot_product_3 = velocity_field[index*3+1];
        double feq_3 = density_field_19 * (dot_product_3 * (c_s_2_inv + dot_product_3 * c_s_4_inv) + norm_square_cs2);
        previous_particle_distributions[scalar_index(x, y, z, 3)] = omtauinv * particle_distributions[scalar_index(x, y, z, 3)] +  feq_3;

        // 7
        double dot_product_4 = -velocity_field[index*3+1];
        double feq_4 = density_field_19 * (dot_product_4 * (c_s_2_inv + dot_product_4 * c_s_4_inv) + norm_square_cs2);
        previous_particle_distributions[scalar_index(x, y, z, 4)] = omtauinv * particle_distributions[scalar_index(x, y, z, 4)] + feq_4;

        // 7
        double dot_product_5 = velocity_field[index*3+2];
        double feq_5 = density_field_19 * (dot_product_5 * (c_s_2_inv + dot_product_5 * c_s_4_inv) + norm_square_cs2);
        previous_particle_distributions[scalar_index(x, y, z, 5)] = omtauinv * particle_distributions[scalar_index(x, y, z, 5)] + feq_5;

        // 7
        double dot_product_6 = -velocity_field[index*3+2];
        double feq_6 = density_field_19 * (dot_product_6 * (c_s_2_inv + dot_product_6 * c_s_4_inv) + norm_square_cs2);
        previous_particle_distributions[scalar_index(x, y, z, 6)] = omtauinv * particle_distributions[scalar_index(x, y, z, 6)] + feq_6;

            // 9
            double dot_product_7 = velocity_field[index*3] + velocity_field[index*3+1] + velocity_field[index*3+2];
            double feq_7 = density_field_172 * (dot_product_7 * (c_s_2_inv + dot_product_7 * c_s_4_inv) + norm_square_cs2);
            previous_particle_distributions[scalar_index(x, y, z, 7)] = omtauinv * particle_distributions[scalar_index(x, y, z, 7)] + feq_7;

            // 9
            double dot_product_8 = -velocity_field[index*3] - velocity_field[index*3+1] - velocity_field[index*3+2];
            double feq_8 = density_field_172 * (dot_product_8 * (c_s_2_inv + dot_product_8 * c_s_4_inv) + norm_square_cs2);
            previous_particle_distributions[scalar_index(x, y, z, 8)] = omtauinv * particle_distributions[scalar_index(x, y, z, 8)] + feq_8;

            // 9
            double dot_product_9 = velocity_field[index*3] + velocity_field[index*3+1] - velocity_field[index*3+2];
            double feq_9 = density_field_172 * (dot_product_9 * (c_s_2_inv + dot_product_9 * c_s_4_inv) + norm_square_cs2);
            previous_particle_distributions[scalar_index(x, y, z, 9)] = omtauinv * particle_distributions[scalar_index(x, y, z, 9)] + feq_9;

            // 9
            double dot_product_10 = -velocity_field[index*3] - velocity_field[index*3+1] + velocity_field[index*3+2];
            double feq_10 = density_field_172 * (dot_product_10 * (c_s_2_inv + dot_product_10 * c_s_4_inv) + norm_square_cs2);
            previous_particle_distributions[scalar_index(x, y, z, 10)] = omtauinv * particle_distributions[scalar_index(x, y, z, 10)] + feq_10;

            // 9
            double dot_product_11 = velocity_field[index*3] - velocity_field[index*3+1] + velocity_field[index*3+2];
            double feq_11 = density_field_172 * (dot_product_11 * (c_s_2_inv + dot_product_11 * c_s_4_inv) + norm_square_cs2);
            previous_particle_distributions[scalar_index(x, y, z, 11)] = omtauinv * particle_distributions[scalar_index(x, y, z, 11)] + feq_11;

            // 9
            double dot_product_12 = -velocity_field[index*3] + velocity_field[index*3+1] - velocity_field[index*3+2];
            double feq_12 = density_field_172 * (dot_product_12 * (c_s_2_inv + dot_product_12 * c_s_4_inv) + norm_square_cs2);
            previous_particle_distributions[scalar_index(x, y, z, 12)] = omtauinv * particle_distributions[scalar_index(x, y, z, 12)] + feq_12;

            // 9
            double dot_product_13 = -velocity_field[index*3] + velocity_field[index*3+1] + velocity_field[index*3+2];
            double feq_13 = density_field_172 * (dot_product_13 * (c_s_2_inv + dot_product_13 * c_s_4_inv) + norm_square_cs2);
            previous_particle_distributions[scalar_index(x, y, z, 13)] = omtauinv * particle_distributions[scalar_index(x, y, z, 13)] + feq_13;

            // 9
        double dot_product_14 = velocity_field[index*3] - velocity_field[index*3+1] - velocity_field[index*3+2];
        double feq_14 = density_field_172 * (dot_product_14 * (c_s_2_inv + dot_product_14 * c_s_4_inv) + norm_square_cs2);
        previous_particle_distributions[scalar_index(x, y, z, 14)] = omtauinv * particle_distributions[scalar_index(x, y, z, 14)] + feq_14;
      }
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

  start_run(stream_profiler);
  stream();
  end_run(stream_profiler);
}

int main(int argc, char const *argv[]) {
  if(system("rm -rf output") == -1){std::cout<<"UNSUCCESSFULLY REMOVED OUTPUT FOLDER"<<std::endl;}
  if(system("mkdir output")== -1){std::cout<<"UNSUCCESSFULLY CREATED OUTPUT FOLDER"<<std::endl;}

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
  return 0;
}
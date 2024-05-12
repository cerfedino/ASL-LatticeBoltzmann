#include <iostream>
#include <stdio.h>
#include "include/utils.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include "main.h"

//for now one can define OUTPUT and TIMING
#define OUTPUT

inline bool file_exists(const std::string &name) {
  struct stat buffer;
  return (stat(name.c_str(), &buffer) == 0);
}

int time_lbm = 0;

inline int scalar_index(int x, int y, int z) {
  return (z * NX * NY) + (y * NX) + x;
}
inline int scalar_index(int x, int y, int z, int w) {
  return (x + y * NX + z * NX * NY + w * NX * NY * NZ);
}

// Lattice directions using D3DQ15. assumed speed of sound c_s = 1/sqrt(3).
vector_3_int *directions;
double *weights;
int *reverse_indexes;
double gamma_dot;
std::string velocity_set = "D3Q15";


void initialise() {
  int box_flatten_length = NX * NY * NZ;
  int distributions_flatten_length = box_flatten_length * direction_size;
  density_field = new double[box_flatten_length];
  velocity_field = new vector_3_double[box_flatten_length];
  previous_particle_distributions = new double[distributions_flatten_length];
  particle_distributions = new double[distributions_flatten_length];
  for (int i = 0; i < NX * NY * NZ; i++) {
    density_field[i] = 1;
  }
  for (int x = 0; x < NX; x++) {
    for (int y = 0; y < NY; y++) {
      for (int z = 0; z < NZ; z++) {
        for (int i = 0; i < direction_size; i++) {
          previous_particle_distributions[scalar_index(x, y, z, i)] =
              weights[i];
          particle_distributions[scalar_index(x, y, z, i)] = weights[i];
        }
      }
    }
  }
}

double calculate_feq(int x, int y, int z, int i) {
  double dot_product =
      (double)velocity_field[scalar_index(x, y, z)].x *
          (double)directions[i].x +
      (double)velocity_field[scalar_index(x, y, z)].y *
          (double)directions[i].y +
      (double)velocity_field[scalar_index(x, y, z)].z * (double)directions[i].z;
  // Equation 3.4 with c_s^2 = 1/3
  double feq =
      weights[i] * density_field[scalar_index(x, y, z)] *
      (1.0 + dot_product / (c_s * c_s) +
       dot_product * dot_product / (2 * c_s * c_s * c_s * c_s) -
       norm_square(velocity_field[scalar_index(x, y, z)]) / (2 * c_s * c_s));
  return feq;
}

double calculate_feq(int x, int y, int z, int i, double u_le_x) {
  double dot_product =
      (velocity_field[scalar_index(x, y, z)].x + u_le_x) * directions[i].x +
      velocity_field[scalar_index(x, y, z)].y * directions[i].y +
      velocity_field[scalar_index(x, y, z)].z * directions[i].z;
  double norm_square =
      (velocity_field[scalar_index(x, y, z)].x + u_le_x) *
          (velocity_field[scalar_index(x, y, z)].x + u_le_x) +
      velocity_field[scalar_index(x, y, z)].y * directions[i].y +
      velocity_field[scalar_index(x, y, z)].z * directions[i].z;
  // Equation 3.4 with c_s^2 = 1/3
  double feq = weights[i] * density_field[scalar_index(x, y, z)] *
               (1.0 + dot_product / (c_s * c_s) +
                dot_product * dot_product / (2 * c_s * c_s * c_s * c_s) -
                norm_square / (2 * c_s * c_s));
  return feq;
}

void set_velocity_set(std::string velocity_set) {
  direction_size = 15;
  directions = new vector_3_int[15]{
      {0, 0, 0},   {1, 0, 0},  {-1, 0, 0},  {0, 1, 0},    {0, -1, 0},
      {0, 0, 1},   {0, 0, -1}, {1, 1, 1},   {-1, -1, -1}, {1, 1, -1},
      {-1, -1, 1}, {1, -1, 1}, {-1, 1, -1}, {-1, 1, 1},   {1, -1, -1}};
  weights = new double[15]{2.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,
                           1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,  1.0 / 72.0,
                           1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0,
                           1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0};
  lookup_reverse();
}

void lookup_reverse() {
  reverse_indexes = new int[direction_size];
  for (int i = 0; i < direction_size; i++) {
    for (int j = 0; j < direction_size; j++) {
      if (directions[i].x == -directions[j].x &&
          directions[i].y == -directions[j].y &&
          directions[i].z == -directions[j].z) {
        reverse_indexes[i] = j;
      }
    }
  }
}

void set_velocity(int x_field, int y_field, int z_field, double u_x, double u_y,
                  double u_z) {
  velocity_field[scalar_index(x_field, y_field, z_field)].x = u_x;
  velocity_field[scalar_index(x_field, y_field, z_field)].y = u_y;
  velocity_field[scalar_index(x_field, y_field, z_field)].z = u_z;
}

void set_density(int x_field, int y_field, int z_field, double density) {
  density_field[scalar_index(x_field, y_field, z_field)] = density;
}

void compute_density_momentum_moment() {
  for (int x = 0; x < NX; x++) {
    for (int y = 0; y < NY; y++) {
      for (int z = 0; z < NZ; z++) {

        double new_density = 0;
        vector_3_double u;
        
        for (int i = 0; i < direction_size; i++) {
          new_density += particle_distributions[scalar_index(x, y, z, i)];
          u.x += particle_distributions[scalar_index(x, y, z, i)] *
                 directions[i].x;
          u.y += particle_distributions[scalar_index(x, y, z, i)] *
                 directions[i].y;
          u.z += particle_distributions[scalar_index(x, y, z, i)] *
                 directions[i].z;
        }

        density_field[scalar_index(x, y, z)] = new_density;
        velocity_field[scalar_index(x, y, z)].x = u.x / new_density;
        velocity_field[scalar_index(x, y, z)].y = u.y / new_density;
        velocity_field[scalar_index(x, y, z)].z = u.z / new_density;
      }
    }
  }
}

void stream() {
  for (int x = 0; x < NX; x++) {
    for (int y = 0; y < NY; y++) {
      for (int z = 0; z < NZ; z++) {
        for (int i = 0; i < direction_size; i++) {
          // Periodic boundary conditions taken from Taylor green in Chapter 13.
            int xmd = (NX + x - (int)directions[i].x) % NX;
            int ymd = (NY + y - (int)directions[i].y) % NY;
            int zmd = (NZ + z - (int)directions[i].z) % NZ;
            particle_distributions[scalar_index(x, y, z, i)] =
                previous_particle_distributions[scalar_index(xmd, ymd, zmd, i)];
          
        }
      }
    }
  }
}

void collision() { // Performs the collision step.
  const double tauinv = 1.0 / tau;
  const double omtauinv = 1.0 - tauinv; // 1 - 1/tau
  for (int x = 0; x < NX; x++) {
    for (int y = 0; y < NY; y++) {
      for (int z = 0; z < NZ; z++) {
        for (int i = 0; i < direction_size; i++) {
          double feq = calculate_feq(x, y, z, i);
          // Equation 3.9
          previous_particle_distributions[scalar_index(x, y, z, i)] =
              omtauinv * particle_distributions[scalar_index(x, y, z, i)] +
              tauinv * feq;
        }
      }
    }
  }
}

void perform_timestep() {
  time_lbm++;
  compute_density_momentum_moment();
  collision();
  stream();
}

int main(int argc, char **argv) {

  int save_every, n_steps;
  int success = setup(save_every, c_s, tau, n_steps, gamma_dot);

  if (success != 1) {
    std::cout << "Setup failed!" << std::endl;
    return 0;
  }

  unsigned long long start_cycle, end_cycle;
  time_t start_sec, end_sec;

  #ifdef TIMING
    asm volatile("RDTSC" : "=A"(start_cycle));
    time(&start_sec);
  #endif
  // setup LBM
  set_velocity_set(velocity_set);
  initialise();

  double viscosity = c_s * c_s * (tau - 0.5);
  
  #ifdef OUTPUT
    output_lbm_data("output/0.csv", true);
    output_indices_file();
  #endif

  int scale = 1;
  int runs = n_steps * scale * scale * scale;

  // start simulation
  for (int i = 0; i < runs; i = i + 1) {
    perform_timestep();
    
    #ifdef OUTPUT
    if ( (i + 1) % save_every == 0) {
      double percentage = (double)(i + 1) / (double)(runs) * 100.0;
      std::cout << "Saving data - " << (i + 1) << "/" << runs << " ("
                << percentage << "%)" << '\n';
      output_lbm_data("output/" + std::to_string(i + 1) + ".csv", true);
    }
    #endif

  }
  std::cout << std::endl;
  free_up();
  
  #ifdef TIMING 
    asm volatile("RDTSC" : "=A"(end_cycle));
    time(&end_sec);
    printf("Cycles taken: %llu (%ld seconds)\n", end_cycle - start_cycle,
           end_sec - start_sec);
  
  #endif
  
  return 0;
}

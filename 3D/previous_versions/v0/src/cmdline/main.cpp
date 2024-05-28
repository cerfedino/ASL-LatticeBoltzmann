#include <iostream>
#include <stdio.h>
// Now Linux only.
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include "profile.h"
#include "utils.h"
#include <iostream>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <fstream>

// the variables are in utils.h
using namespace std;

// for now one can define OUTPUT, TIMING
#define OUTPUT

profiler *compute_density_momentum_profiler, *collision_profiler;

int time_lbm = 0;


int NX, NY, NZ, NT;


double *density_field;
vector_3_double *velocity_field;
double *previous_particle_distributions;
double *particle_distributions;
int direction_size = 15;

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
  velocity_field = (vector_3_double *)malloc(box_flatten_length * sizeof(vector_3_double));

  previous_particle_distributions = (double *)malloc(distributions_flatten_length * sizeof(double));
  particle_distributions = (double *)malloc(distributions_flatten_length * sizeof(double));

  for (int i = 0; i < NX * NY * NZ; i++) {
    density_field[i] = 1;
  }
  for (int x = 0; x < NX; x++) {
    for (int y = 0; y < NY; y++) {
      for (int z = 0; z < NZ; z++) {
        for (int i = 0; i < direction_size; i++) {
          previous_particle_distributions[scalar_index(x, y, z, i)] = weights[i];
          particle_distributions[scalar_index(x, y, z, i)] = weights[i];
        }
      }
    }
  }

  compute_density_momentum_profiler = init_profiler(NX * NY * NZ * direction_size * 7 + NX * NY * NZ * 3, 8 * NX * NY * NZ * direction_size * 7 + (NX * NY * NZ * 3 + 15) * 8);
  collision_profiler = init_profiler((NX * NY * NZ * direction_size * 29 + 2), 8 * (4 * NX * NY * NZ + 15 + 15));
}

double calculate_feq(int x, int y, int z, int i) {
  double dot_product = (double)velocity_field[scalar_index(x, y, z)].x * (double)directions[i].x + (double)velocity_field[scalar_index(x, y, z)].y * (double)directions[i].y + (double)velocity_field[scalar_index(x, y, z)].z * (double)directions[i].z;

  double feq = weights[i] * density_field[scalar_index(x, y, z)] * (1.0 + dot_product / (c_s * c_s) + dot_product * dot_product / (2 * c_s * c_s * c_s * c_s) - norm_square(velocity_field[scalar_index(x, y, z)]) / (2 * c_s * c_s));
  return feq;
}

double calculate_feq(int x, int y, int z, int i, double u_le_x) {
  double dot_product = (velocity_field[scalar_index(x, y, z)].x + u_le_x) * directions[i].x + velocity_field[scalar_index(x, y, z)].y * directions[i].y + velocity_field[scalar_index(x, y, z)].z * directions[i].z;
  double norm_square = (velocity_field[scalar_index(x, y, z)].x + u_le_x) * (velocity_field[scalar_index(x, y, z)].x + u_le_x) + velocity_field[scalar_index(x, y, z)].y * directions[i].y + velocity_field[scalar_index(x, y, z)].z * directions[i].z;
  // Equation 3.4 with c_s^2 = 1/3
  double feq = weights[i] * density_field[scalar_index(x, y, z)] * (1.0 + dot_product / (c_s * c_s) + dot_product * dot_product / (2 * c_s * c_s * c_s * c_s) - norm_square / (2 * c_s * c_s));
  return feq;
}

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

  reverse_indexes = (int *)malloc(direction_size * sizeof(int));
  for (int i = 0; i < direction_size; i++) {
    for (int j = 0; j < direction_size; j++) {
      if (directions[i].x == -directions[j].x && directions[i].y == -directions[j].y && directions[i].z == -directions[j].z) {
        reverse_indexes[i] = j;
      }
    }
  }
}

void set_velocity(int x_field, int y_field, int z_field, double u_x, double u_y, double u_z) {
  velocity_field[scalar_index(x_field, y_field, z_field)].x = u_x;
  velocity_field[scalar_index(x_field, y_field, z_field)].y = u_y;
  velocity_field[scalar_index(x_field, y_field, z_field)].z = u_z;
}

void set_density(int x_field, int y_field, int z_field, double density) { density_field[scalar_index(x_field, y_field, z_field)] = density; }

void compute_density_momentum_moment() {
	int scalar_index_curr = 0;
	//we got some terrible memory accesses here... 
	// we gotta do some blocking most probably
	// oof
  for (int x = 0; x < NX; x++) {
    for (int y = 0; y < NY; y++) {
      for (int z = 0; z < NZ; z++) {

        double new_density = 0, new_density_1 = 0, new_density_2 = 0, new_density_3 = 0, new_density_4 = 0, new_density_5 = 0;
        //vector_3_double u, u1, u2, u3, u4, u5;
		double x_1=0, x_2=0, x_3=0, x_4=0, x_5=0;
		double y_1=0, y_2=0, y_3=0, y_4=0, y_5=0;
		double z_1=0, z_2=0, z_3=0, z_4=0, z_5=0;

        for (int i = 0; i < direction_size; i++) {
		  scalar_index_curr = scalar_index(x, y, z, i);
          new_density_1 += particle_distributions[scalar_index_curr];
          x_1 += particle_distributions[scalar_index_curr] * directions[i].x;
          y_1 += particle_distributions[scalar_index_curr] * directions[i].y;
          z_1 += particle_distributions[scalar_index_curr] * directions[i].z;
		  scalar_index_curr+=NX * NY * NZ;
		  i++;

		  new_density_2 += particle_distributions[scalar_index_curr];
          x_2 += particle_distributions[scalar_index_curr] * directions[i].x;
          y_2 += particle_distributions[scalar_index_curr] * directions[i].y;
          z_2 += particle_distributions[scalar_index_curr] * directions[i].z;
		  scalar_index_curr+=NX * NY * NZ;
		  i++;

		  new_density_3 += particle_distributions[scalar_index_curr];
          x_3 += particle_distributions[scalar_index_curr] * directions[i].x;
          y_3 += particle_distributions[scalar_index_curr] * directions[i].y;
          z_3 += particle_distributions[scalar_index_curr] * directions[i].z;
		  scalar_index_curr+=NX * NY * NZ;
		  i++;

		  new_density_4 += particle_distributions[scalar_index_curr];
          x_4 += particle_distributions[scalar_index_curr] * directions[i].x;
          y_4 += particle_distributions[scalar_index_curr] * directions[i].y;
          z_4 += particle_distributions[scalar_index_curr] * directions[i].z;
		  scalar_index_curr+=NX * NY * NZ;
		  i++;

		  new_density_5 += particle_distributions[scalar_index_curr];
          x_5 += particle_distributions[scalar_index_curr] * directions[i].x;
          y_5 += particle_distributions[scalar_index_curr] * directions[i].y;
          z_5 += particle_distributions[scalar_index_curr] * directions[i].z;
		
		}

		new_density = new_density_1 + new_density_2 + new_density_3+new_density_4+new_density_5;
		double x_sum = x_1+x_2+x_3+x_4+x_5;
		double y_sum = y_1+y_2+y_3+y_4+y_5;
		double z_sum = z_1+z_2+z_3+z_4+z_5;

        density_field[scalar_index(x, y, z)] = new_density;
        velocity_field[scalar_index(x, y, z)].x = x_sum / new_density;
        velocity_field[scalar_index(x, y, z)].y = y_sum / new_density;
        velocity_field[scalar_index(x, y, z)].z = z_sum / new_density;
      }
    }
  }

}

void stream() {
    for(int x = 0; x < NX; x++) {
        for(int y = 0; y < NY; y++) {
            for(int z = 0; z < NZ; z++) {
                for(int i = 0; i < direction_size; i++) {
                    //Periodic boundary conditions taken from Taylor green in Chapter 13.
                    if (boundary_conditions == "periodic") {
                        //X position Minus the Direction (xmd) applies to y and z.
                        int xmd = (NX + x - (int)directions[i].x) % NX;
                        int ymd = (NY + y - (int)directions[i].y) % NY;
                        int zmd = (NZ + z - (int)directions[i].z) % NZ;
                        particle_distributions[scalar_index(x,y,z,i)] = previous_particle_distributions[scalar_index(xmd,ymd,zmd,i)];
                        //Equation 3.10 with periodic boundary conditions.
                    } else if(boundary_conditions == "couette") {

                        //Using Bounce back approach for top & bottom walls..
                        if (y==0 && directions[i].y == 1) {
                            //Bottom Wall.
                            //Equation 5.27 from LBM Principles and Practice.
                            particle_distributions[scalar_index(x,y,z,i)]=previous_particle_distributions[scalar_index(x,y,z,reverse_indexes[i])];
                        } else if (y==NY-1 && directions[i].y == -1) {
                            //Top wall
                            //Equation 5.28 from LBM Principles and Practice.
                            //coefficients of Equation 5.28 calculated from footnote 17.
                            double u_max = 0.1;
                            particle_distributions[scalar_index(x,y,z,i)]=previous_particle_distributions[scalar_index(x,y,z,reverse_indexes[i])]+directions[i].x * 2 * weights[i] / (c_s * c_s) * u_max;
                            //particle_distributions[scalar_index(x,y,z,(4-1))]=previous_particle_distributions[scalar_index(x,y,z,(2-1))];
                            //particle_distributions[scalar_index(x,y,z,(7-1))]=previous_particle_distributions[scalar_index(x,y,z,(5-1))]-(1.0/6.0)*u_max;
                            //particle_distributions[scalar_index(x,y,z,(8-1))]=previous_particle_distributions[scalar_index(x,y,z,(6-1))]+(1.0/6.0)*u_max;
                        } else {
                            //Chapter 13 Taylor-Green periodicity from same book.
                            int xmd = (NX + x - directions[i].x) % NX;

                            //int ymd = (NY + y - (int)directions[i].y) % NY;
                            int ymd = y - directions[i].y;
                            int zmd = (NZ + z - directions[i].z) % NZ;
                            particle_distributions[scalar_index(x,y,z,i)] = previous_particle_distributions[scalar_index(xmd,ymd,zmd,i)];
                        }
                    } else if(boundary_conditions == "lees_edwards") {
                            double d_x = gamma_dot * (double)NY * (double)time_lbm;
                            int d_x_I = d_x;
                            double d_x_R = d_x - d_x_I;
                            d_x_I = d_x_I % NX;
                            double s_1 = d_x_R;
                            double s_2 = 1.0 - d_x_R;
                            //s_1=0;s_2=1;
                            int ymd = (NY + y - directions[i].y) % NY;
                            int zmd = (NZ + z - directions[i].z) % NZ;
                            if(y==0 && directions[i].y == 1) {
                                //Bottom Wall.
                                //Equation (17) from Less-Edwards boundary conditions for lattice Boltzmann suspension simulations
                                //by Eric Lorenz and Alfons G. Hoekstra
                                    int x_pos;
                                    x_pos = (x+d_x_I + 2 * NX - directions[i].x) % NX;
                                    int x_shifted;
                                    x_shifted = (x+d_x_I + 2 * NX +1 - directions[i].x) % NX;
                                    double galilean_transformation_pos = previous_particle_distributions[scalar_index(x_pos,ymd,zmd,i)] + calculate_feq(x_pos,ymd,zmd,i,-1 * gamma_dot * (double)NY) - calculate_feq(x_pos,ymd,zmd,i,0);
                                    double galilean_transformation_shift = previous_particle_distributions[scalar_index(x_shifted,ymd,zmd,i)] + calculate_feq(x_shifted,ymd,zmd,i,-1 * gamma_dot * (double)NY) - calculate_feq(x_shifted,ymd,zmd,i,0);
                                    //Equation (18) from same paper.
                                    particle_distributions[scalar_index(x,y,z,i)]= s_1 * galilean_transformation_shift +
                                                                                   s_2 * galilean_transformation_pos;
                                } else if(y==NY-1 && directions[i].y == -1) {
                                    //Equation (17) from Less-Edwards boundary conditions for lattice Boltzmann suspension simulations
                                    //by Eric Lorenz and Alfons G. Hoekstra
                                    //Top Wall.
                                    int x_pos;
                                    x_pos = (x-d_x_I + 2 * NX - directions[i].x) % NX;
                                    int x_shifted;
                                    x_shifted = (x-d_x_I + 2 * NX -1 - directions[i].x) % NX;
                                    double galilean_transformation_pos = previous_particle_distributions[scalar_index(x_pos,ymd,zmd,i)] + calculate_feq(x_pos,ymd,zmd,i,1 * gamma_dot * (double)NY) - calculate_feq(x_pos,ymd,zmd,i,0);
                                    double galilean_transformation_shift = previous_particle_distributions[scalar_index(x_shifted,ymd,zmd,i)] + calculate_feq(x_shifted,ymd,zmd,i,1 * gamma_dot * (double)NY) - calculate_feq(x_shifted,ymd,zmd,i,0);
                                    //Equation (18) from same paper.
                                    particle_distributions[scalar_index(x,y,z,i)]= s_1 * galilean_transformation_shift +
                                                                                   s_2 * galilean_transformation_pos;
                                } else {
                                    //Interior points. (Enforce periodicity along x axis);.
                                    int xmd = (NX + x - directions[i].x) % NX;
                                    particle_distributions[scalar_index(x,y,z,i)] = previous_particle_distributions[scalar_index(xmd,ymd,zmd,i)];
                                }


                    }
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
          previous_particle_distributions[scalar_index(x, y, z, i)] = omtauinv * particle_distributions[scalar_index(x, y, z, i)] + tauinv * feq;
        }
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

  stream();
}

int main(int argc, char const *argv[]) {
  if (argc == 5) {
    NX = atoi(argv[1]);
    NY = atoi(argv[2]);
    NZ = atoi(argv[3]);
    NT = atoi(argv[4]);
  } else if (argc != 1) {
    printf("Usage: %s [NX NY NZ NT]\n", argv[0]);
    return 1;
  }
  
  system("rm -rf output");
  system("mkdir output");

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
  //std::cout<<viscosity<<std::endl;

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
      if (((i + 1) % save_every) == 0) {
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
#endif
  return 0;
}
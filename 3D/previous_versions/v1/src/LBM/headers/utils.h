#ifndef UTILS_H
#define UTILS_H

#define save_every 20
#define c_s 0.577350269
#define velocity_set "D3Q15"
#define nu 0.16666666666
#define tau 0.6
#define gamma_dot 0.01
#define boundary_conditions "couette"

#include <string>

struct vector_3_int {
  int x = 0;
  int y = 0;
  int z = 0;
  vector_3_int(int x_val, int y_val, int z_val) {
    x = x_val;
    y = y_val;
    z = z_val;
  }
};

struct vector_3_double {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

inline double norm_square(vector_3_double d) { return (d.x * d.x + d.y * d.y + d.z * d.z); }

extern double *density_field;
extern vector_3_double *velocity_field;
extern double *previous_particle_distributions;
extern double *particle_distributions;
extern int direction_size;
extern int NX, NY, NZ, NT;

inline void free_up() {
  free(density_field);
  free(particle_distributions);
}

void output_array(double *array);
void output_density();
void output_velocity();
void output_indices_file();
void output_lbm_data(std::string filename, bool header);
void output_f_array(double *f_array, int z_index);
#endif
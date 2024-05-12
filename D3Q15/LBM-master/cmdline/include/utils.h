#ifndef UTILS_H
#define UTILS_H

#include <string>


struct vector_3_int{
  int x=0;
  int y=0; 
  int z=0;
vector_3_int(int x_val, int y_val, int z_val)  {x=x_val;y=y_val;z=z_val;}
};

struct vector_3_double{
  double x=0.0;
  double y=0.0; 
  double z=0.0;
};

inline double norm_square(vector_3_double d) {
	return (d.x * d.x + d.y * d.y + d.z * d.z);
}


extern int NX;
extern int NY;
extern int NZ;

extern double c_s;
extern double nu;
extern double tau;

extern double *density_field;
extern vector_3_double *velocity_field;
extern double *previous_particle_distributions;
extern double *particle_distributions;
extern int direction_size;

int setup(int& save_every, double& c_s, double& tau, int& n_steps, double& gamma_dot);

inline void free_up() {
  free(density_field);
  free(particle_distributions);
}

void output_array(double *array);
void output_density();
void output_velocity();
void output_indices_file();
void output_lbm_data(std::string filename, bool header);
void output_f_array(double *f_array, int z_index) ;
#endif
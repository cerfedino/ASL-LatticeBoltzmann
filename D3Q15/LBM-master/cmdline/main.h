void set_velocity(int x_field, int y_field, int z_field, double u_x, double u_y,
                  double u_z); // Set velocity at position in velocity field.
void set_density(int x_field, int y_field, int z_field,
                 double density); // Set density at position in density field.
double calculate_feq(int i, int j, int k, int w);
double calculate_feq(int i, int j, int k, int w, double u_le_x);

void output_velocity();
void output_test();
void output_indices_file();
void compute_density_momentum_moment();
void stream(); // Stream the current equilibrium distribution to the next
               // distribution.
void collision(); // Perform the collision step. Assumes delta t / tau = 1.
void perform_timestep(); // Delta t = 1 lattice unit.
void output_lbm_data(std::string filename, bool header = true);
void lookup_reverse();
void output_f_array(double *f_array, int z_index);
int get_time();
#include "utils.h"
#include <fstream>
#include <iostream>
// Now Linux only.
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

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

int scalar_index(int x, int y, int z) { return (z * NX * NY) + (y * NX) + x; }
int scalar_index(int x, int y, int z, int w) { return (x + y * NX + z * NX * NY + w * NX * NY * NZ); }

void output_array(double *array) {
  std::cout << "x,y,z value" << std::endl;
  for (int x = 0; x < NX; x++) {
    for (int y = 0; y < NY; y++) {
      for (int z = 0; z < NZ; z++) {
        std::cout << x << "," << y << "," << z << ": " << array[scalar_index(x, y, z)] << std::endl;
      }
    }
  }
}

void output_density() { output_array(density_field); }

void output_velocity() {
  int z_index = 0;
  for (int x = 0; x < NX; x++) {
    for (int y = 0; y < NY; y++) {
      std::cout << velocity_field[scalar_index(x, y, z_index)*3] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void output_indices_file() {
  std::ofstream output("output/indices.csv");
  output << "x,y,z" << '\n';
  for (int i = 0; i < NX; i++) {
    for (int j = 0; j < NY; j++) {
      for (int k = 0; k < NZ; k++) {
        output << i << "," << j << "," << k << '\n';
      }
    }
  }
  std::cout << std::endl;
  output.close();
}

void output_f_array(double *f_array, int z_index) {
  for (int i = 0; i < direction_size; i++) {
    std::cout << "ans(:,:," << z_index << "," << (i + 1) << ")" << '\n';
    std::cout << '\n';
    for (int x = 0; x < NX; x++) {
      for (int y = 0; y < NY; y++) {
        std::cout << f_array[scalar_index(x, y, z_index, i)] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << '\n';
  }
  std::cout << std::endl;
}

void output_lbm_data(std::string filename, bool header) {
  std::ofstream output_stream;
  output_stream.open(filename, std::ofstream::out | std::ofstream::app);
  if (header) {
    output_stream << "p,u_x,u_y,u_z" << '\n';
  }
  for (int x = 0; x < NX; x++) {
    for (int y = 0; y < NY; y++) {
      for (int z = 0; z < NZ; z++) {
        output_stream << density_field[scalar_index(x, y, z)] << "," << velocity_field[scalar_index(x, y, z)*3] << "," << velocity_field[scalar_index(x, y, z)*3+1] << "," << velocity_field[scalar_index(x, y, z)*3+2] << '\n';
      }
    }
  }
  output_stream.close();
}

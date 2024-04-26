#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <cmath>
#include <exception>
#include <stdexcept>
#include "include/npy.hpp"
#include <vector>
#include <string>

#define PRINT 1

#define debug_printf(fmt, ...)           \
  do {                                   \
    if (PRINT) {                         \
      fprintf(stdout, fmt, __VA_ARGS__); \
    }                                    \
  } while (0)

#define debug_print(fmt)    \
  do {                      \
    if (PRINT) {            \
      fprintf(stdout, fmt); \
    }                       \
  } while (0)

using namespace std;

const int Nx = 400;    // resolution in x
const int Ny = 100;    // resolution in y
const int rho0 = 100;  // average density
const int tau = 0.6;   // collision timescale
const int Nt = 100;   // number of timesteps

// Lattice speeds / weights
const int NL = 9;
const int idx[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
const int cxs[9] = {0, 0, 1, 1, 1, 0, -1, -1, -1};
const int cys[9] = {0, 1, 1, 0, -1, -1, -1, 0, 1};
const double w[9] = {4 / 9, 1 / 9,  1 / 36, 1 / 9, 1 / 36,
                     1 / 9, 1 / 36, 1 / 9,  1 / 36};  // sums to 1

void meshgrid(int **x_coords, int **y_coords) {
  // Fill coordinate matrices
  for (int i = 0; i < Ny; ++i) {
    for (int j = 0; j < Nx; ++j) {
      x_coords[i][j] = j;
      y_coords[i][j] = i;
    }
  }
}

void roll_array(int i, int cx, int cy, double ***F) {
  // Temporary array to hold rolled values
  double temp[Ny][Nx];

  // Roll operation along axis 1 (cx)
  for (int j = 0; j < Ny; ++j) {
    for (int k = 0; k < Nx; ++k) {
      temp[j][(k + cx + Nx) % Nx] = F[i][j][k];
    }
  }

  // Roll operation along axis 0 (cy)
  for (int j = 0; j < Ny; ++j) {
    for (int k = 0; k < Nx; ++k) {
      F[i][(j + cy + Ny) % Ny][k] = temp[j][k];
    }
  }
}


// TODO if fancy make generic functions for double, float, int and so on

double ***malloc_3d(int x, int y, int z) {
  double ***array = (double ***)malloc(x * sizeof(double **));
  for (int i = 0; i < x; ++i) {
    array[i] = (double **)malloc(y * sizeof(double *));
    for (int j = 0; j < y; ++j) {
      array[i][j] = (double *)malloc(z * sizeof(double));
    }
  }
  return array;
}

int **malloc_2d(int x, int y) {
  int **array = (int **)malloc(x * sizeof(int *));
  for (int i = 0; i < x; ++i) {
    array[i] = (int *)malloc(y * sizeof(int));
  }
  return array;
}

double **malloc_2d_double(int x, int y) {
  double **array = (double **)malloc(x * sizeof(double *));
  for (int i = 0; i < x; ++i) {
    array[i] = (double *)malloc(y * sizeof(double));
  }
  return array;
}

void save_npy_3d_double(double ***array, int x, int y, int z, string filename) {
  
  npy::npy_data_ptr<double> d;
  d.data_ptr = &array[0][0][0];
  d.shape = {(unsigned long)x, (unsigned long)y, (unsigned long)z};
  d.fortran_order = false; // optional

  const std::string path{filename};
  npy::write_npy(path, d);
}

void save_npy_2d_double(double **array, int x, int y, string filename) {
  npy::npy_data_ptr<double> d;
  d.data_ptr = &array[0][0];
  d.shape = {(unsigned long)x, (unsigned long)y};
  d.fortran_order = false; // optional

  const std::string path{filename};
  npy::write_npy(path, d);
}

void save_npy_2d_int(int **array, int x, int y, string filename) {
  debug_printf("Saving %s\n", filename.c_str());

  debug_printf("Pointer: %p\n", &array[0][0]);

  // print first 10 values
  for (int i = 0; i < 10; i++) {
    debug_printf("%d ", array[0][i]);
  }

  debug_print("\n");

  npy::npy_data_ptr<int> d;
  d.data_ptr = &array[0][0];
  d.shape = {(unsigned long)x, (unsigned long)y};
  d.fortran_order = false; // optional

  debug_printf("Creating %s\n", filename.c_str());

  const std::string path{filename};
  npy::write_npy(path, d);

  debug_printf("Saved %s\n", filename.c_str());
}

string prepare_output() {
  // check if output folder exists in the current folder
  // if not create it
  if (system("ls | grep output") != 0) {
    debug_print("Creating main output folder\n");
    system("mkdir output");

    // check if output folder was created
    if (system("ls | grep output") != 0) {
      debug_print("Could not create output folder\n");
      //throw runtime_error("Could not create output folder");
    }
  }

  // create new folder of format YYYY_MM_DD_HH_MM_SS
  time_t now = time(0);
  tm *ltm = localtime(&now);

  char folder_name[100];
  sprintf(folder_name, "output/%d_%d_%d_%d_%d_%d", 1900 + ltm->tm_year,
          1 + ltm->tm_mon, ltm->tm_mday, ltm->tm_hour, ltm->tm_min,
          ltm->tm_sec);

  char command[100];
  sprintf(command, "mkdir output/%d_%d_%d_%d_%d_%d", 1900 + ltm->tm_year,
          1 + ltm->tm_mon, ltm->tm_mday, ltm->tm_hour, ltm->tm_min,
          ltm->tm_sec);

  // create folder in output folder
  //system("cd output");

  // "mkdir" + folder_name

  system((const char *)command);
  //system("cd ..");

  if (system(folder_name) != 0) {
    debug_printf("Could not create output folder with name %s\n", folder_name);
    /*throw runtime_error("Could not create output folder with name " +
                        string(folder_name));*/
  }

  return string(folder_name);
}

// TODO CHECK TYPES IF ALL CORRECT
// TODO arrays need to be initialized likely with 0 values
int main() {
  string folder_name = prepare_output(); // TODO delete empty folders

  debug_printf("Output folder: %s\n", folder_name.c_str());

  // Lattice Boltzmann Simulation in 2D
  debug_print("Starting\n");

  // Initial conditions
  // distribution function

  // TODO creating array directly segfaults because of size
  double ***F = malloc_3d(Nx, Ny, NL);
  debug_print("Initializing\n");

  // srand(42);  // some seed

  debug_print("Seed set\n");

  // Initialize F
  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < NL; k++) {
        F[i][j][k] = 1 + 0.01 * 0.5;  // TODO randomize
      }
    }
  }

  save_npy_3d_double(F, Nx, Ny, NL, folder_name + "/F.npy");

  // TODO maybe also malloc at some point
  int **x_coords = malloc_2d(Nx, Ny);
  int **y_coords = malloc_2d(Nx, Ny);

  meshgrid(x_coords, y_coords);

  // F[:,:,3] += 2 * (1+0.2*np.cos(2*np.pi*X/Nx*4))
  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
      F[i][j][3] += 2 * (1 + 0.2 * cos(2 * M_PI * x_coords[i][j] / Nx * 4));
    }
  }

  // rho = np.sum(F, axis=2)
  int rho = 0;
  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < NL; k++) {
        rho += F[i][j][k];
      }
    }
  }

  // 	for i in idxs:		F[:,:,i] *= rho0 / rho
  for (int i = 0; i < NL; i++) {
    for (int j = 0; j < Nx; j++) {
      for (int k = 0; k < Ny; k++) {
        F[j][k][i] *= rho0 / rho;
      }
    }
  }

  // Cylinder boundary
  int **cylinderX = malloc_2d(Nx, Ny);
  int **cylinderY = malloc_2d(Nx, Ny);
  int **cylinder = malloc_2d(Nx, Ny);

  meshgrid(cylinderX, cylinderY);

  // cylinder = (X - Nx/4)**2 + (Y - Ny/2)**2 < (Ny/4)**2
  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
      if ((pow(x_coords[i][j] - Nx / 4, 2) + pow(y_coords[i][j] - Ny / 2, 2)) <
          pow(Ny / 4, 2)) {
        cylinder[i][j] = 1;
      } else {
        cylinder[i][j] = 0;
      }
    }
  }

  save_npy_2d_int(cylinder, Nx, Ny, folder_name + "/cylinder.npy");

  debug_print("Starting simulation\n");

  // Simulation
  for (int i = 0; i < Nt; i++) {
    printf("Timestep %d\n", i);

    // Drift IDK about this one
    for (int j = 0; j < NL; j++) {
      roll_array(j, cxs[j], cys[j], F);
    }

    // bndryF = F[cylinder,:]
    // its 2d of size 1941x9 but no idea how this is calculated ??????
    double ***bndryF = malloc_3d(Nx, Ny, NL);
    for (int j = 0; j < Nx; j++) {
      for (int k = 0; k < Ny; k++) {
        if (cylinder[j][k] == 1) {
          for (int l = 0; l < NL; l++) {
            bndryF[j][k][l] = F[j][k][l];
          }
        }
      }
    }

    // translate bndryF = bndryF[:,[0,5,6,7,8,1,2,3,4]]
    for (int j = 0; j < Nx; j++) {
      for (int k = 0; k < Ny; k++) {
        // TODO bnryF[:,[0,5,6,7,8,1,2,3,4]]
      }
    }

    // np.sum(F,2) sum over 2nd axis
    double **rho = malloc_2d_double(Nx, Ny);
    for (int j = 0; j < Nx; j++) {
      for (int k = 0; k < Ny; k++) {
        for (int l = 0; l < NL; l++) {
          rho[j][k] += F[j][k][l];
        }
      }
    }

    // ux = np.sum(F * cxs, 2) / rho
    double **ux = malloc_2d_double(Nx, Ny);
    for (int j = 0; j < Nx; j++) {
      for (int k = 0; k < Ny; k++) {
        for (int l = 0; l < NL; l++) {
          ux[j][k] += F[j][k][l] * cxs[l];
        }
        ux[j][k] /= rho[j][k];
      }
    }

    // uy = np.sum(F * cys, 2) / rho
    double **uy = malloc_2d_double(Nx, Ny);
    for (int j = 0; j < Nx; j++) {
      for (int k = 0; k < Ny; k++) {
        for (int l = 0; l < NL; l++) {
          uy[j][k] += F[j][k][l] * cys[l];
        }
        uy[j][k] /= rho[j][k];
      }
    }

    double ***Feq = malloc_3d(Nx, Ny, NL);
    // set to zero
    for (int j = 0; j < Nx; j++) {
      for (int k = 0; k < Ny; k++) {
        for (int l = 0; l < NL; l++) {
          Feq[j][k][l] = 0;
        }
      }
    }

    debug_printf("Feq shape: %d %d %d\n", Nx, Ny, NL);

    // for i, cx, cy, w in zip(idxs, cxs, cys, weights):
    //         Feq[:,:,i] = rho * w * ( 1 + 3*(cx*ux+cy*uy)  +
    //         9*(cx*ux+cy*uy)**2/2 - 3*(ux**2+uy**2)/2 )
    for (int k = 0; k < NL; k++) {
      for (int j = 0; j < Nx; j++) {
        for (int l = 0; l < Ny; l++) {
          Feq[j][l][k] =
              rho[j][l] * w[k] *
              (1 + 3 * (cxs[k] * ux[j][l] + cys[k] * uy[j][l]) +
               9 * pow(cxs[k] * ux[j][l] + cys[k] * uy[j][l], 2) / 2 -
               3 * (pow(ux[j][l], 2) + pow(uy[j][l], 2)) / 2);
        }
      }
    }

    // F += -(1.0/tau) * (F - Feq)
    for (int j = 0; j < Nx; j++) {
      for (int k = 0; k < Ny; k++) {
        for (int l = 0; l < NL; l++) {
          F[j][k][l] += -(1.0 / tau) * (F[j][k][l] - Feq[j][k][l]);
        }
      }
    }

    // Apply boundary
    // F[cylinder,:] = bndryF
    for (int j = 0; j < Nx; j++) {
      for (int k = 0; k < Ny; k++) {
        if (cylinder[j][k] == 1) {
          for (int l = 0; l < NL; l++) {
            F[j][k][l] = bndryF[j][k][l];
          }
        }
      }
    }

    // vorticity = (np.roll(ux, -1, axis=0) - np.roll(ux, 1, axis=0)) -
    // (np.roll(uy, -1, axis=1) - np.roll(uy, 1, axis=1)) vorticity[cylinder]
    // = np.nan vorticity = np.ma.array(vorticity, mask=cylinder)
    double **vorticity = malloc_2d_double(Nx, Ny);
    for (int j = 0; j < Nx; j++) {
      for (int k = 0; k < Ny; k++) {
        //vorticity[j][k] = (ux[(j + 1 + Nx) % Nx][k] - ux[(j - 1 + Nx) % Nx][k]) -
        //                  (uy[j][(k + 1 + Ny) % Ny] - uy[j][(k - 1 + Ny) % Ny]);
      }
    }

    // Save vorticity
    save_npy_2d_double(vorticity, Nx, Ny, folder_name + "/step_" + to_string(i) + "_vorticity.npy");
  }

  return 0;
}

// to run and compile run: clear && g++ -o 2D/2D 2D/main.cpp -lm && 2D/2D
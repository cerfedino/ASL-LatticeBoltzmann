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
const double rho0 = 100;  // average density
const double tau = 0.6;   // collision timescale
const int Nt = 5000;   // number of timesteps

// Lattice speeds / weights
const int NL = 9;
const double idx[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
const double cxs[9] = {0.0, 0.0, 1.0, 1.0, 1.0, 0.0, -1.0, -1.0, -1.0};
const double cys[9] = {0.0, 1.0, 1.0, 0.0, -1.0, -1.0, -1.0, 0.0, 1.0};
const double weights[9] = {4.0 / 9, 1.0 / 9,  1.0 / 36, 1.0 / 9, 1.0 / 36,
                     1.0 / 9, 1.0 / 36, 1.0 / 9,  1.0 / 36};  // sums to 1

void meshgrid(int **x_coords, int **y_coords) {
  for (int i = 0; i < Ny; ++i) {
    for (int j = 0; j < Nx; ++j) {
      x_coords[i][j] = j;
      y_coords[i][j] = i;
    }
  }
}

void meshgrid(double **x_coords, double **y_coords) {
  for (int i = 0; i < Ny; ++i) {
    for (int j = 0; j < Nx; ++j) {
      x_coords[i][j] = j;
      y_coords[i][j] = i;
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
  
  // convert array to vector  
  // this is needed because npy::write_npy expects a vector and cant deal with pointers to arrays somehow

  vector<double> vec;

  for (int i = 0; i < x; i++) {
    for (int j = 0; j < y; j++) {
      for (int k = 0; k < z; k++) {
        vec.push_back(array[i][j][k]);
      }
    }
  }


  npy::npy_data<double> d;
  d.data = vec;
  d.shape = {(unsigned long)x, (unsigned long)y, (unsigned long)z};
  d.fortran_order = false; // optional

  const std::string path{filename};
  npy::write_npy(path, d);
}

void save_npy_2d_double(double **array, int x, int y, string filename) {

  // convert array to vector  
  // this is needed because npy::write_npy expects a vector and cant deal with pointers to arrays somehow

  vector<double> vec;

  for (int i = 0; i < x; i++) {
    for (int j = 0; j < y; j++) {
      vec.push_back(array[i][j]);
    }
  }


  npy::npy_data<double> d;
  d.data = vec;
  d.shape = {(unsigned long)x, (unsigned long)y};
  d.fortran_order = false; // optional

  const std::string path{filename};
  npy::write_npy(path, d);
}

void save_npy_2d_int(int **array, int x, int y, string filename) {

  // convert array to vector
  // this is needed because npy::write_npy expects a vector and cant deal with pointers to arrays somehow
  // NOTE THIS FUNCTION SHOULD LIKELY NOT BE USED AS INT STORING SOMEHOW DOESNT WORK AND SAVES INCORRECTLY THE DATA

  vector<int> vec;
  for (int i = 0; i < x; i++) {
    for (int j = 0; j < y; j++) {
      vec.push_back(array[i][j]);
    }
  }

  npy::npy_data<int> d;
  d.data = vec;
  d.shape = {(unsigned long)x, (unsigned long)y};
  d.fortran_order = false; // optional

  const std::string path{filename};
  npy::write_npy(path, d);
}


//orig [1 2 3 4 5]
//roll [4 5 1 2 3]
// if called on np.roll(array, 2)
void roll1D(double *array, int size, int shift) {
  double *temp = (double *)malloc(size * sizeof(double));

  for (int i = 0; i < size; i++) {
    temp[(i + shift + size) % size] = array[i];
  }

  for (int i = 0; i < size; i++) {
    array[i] = temp[i];
  }

  free(temp);
}


/*orig
 [[1 2 3]
 [4 5 6]
 [7 8 9]]
roll
 [[3 1 2]
 [6 4 5]
 [9 7 8]]
 if called on np.roll(array, 1, axis=0)
 */
// similar to roll1D
void roll2D(double **array, int height, int width, int shift, int axis){
  double **temp = malloc_2d_double(height, width);

  // matrix is array[Y][X]

  // axis = 0 -> roll along x
  // axis = 1 -> roll along y

  if(axis == 0){
    for (int i = 0; i < height; i++) {
      for (int j = 0; j < width; j++) {
        temp[i][(j + shift + width) % width] = array[i][j];
      }
    }
  } else {
    for (int i = 0; i < height; i++) {
      for (int j = 0; j < width; j++) {
        temp[(i + shift + height) % height][j] = array[i][j];
      }
    }
  }
     

  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      array[i][j] = temp[i][j];
    }
  }

  free(temp);
}

/// @brief This function prepares the output main output folder and creates a new folder with the current date and time
/// @return 
string prepare_output() {
  if (system("ls | grep output") != 0) {
    debug_print("Creating main output folder\n");
    system("mkdir output");

    if (system("ls | grep output") != 0) {
      debug_print("Could not create output folder\n");
    }
  }

  // create new folder of format YYYY_MM_DD_HH_MM_SS
  time_t now = time(0);
  tm *ltm = localtime(&now);

  // its this ugly because i cant be bothered to do const char shenanigans
  char folder_name[100];
  sprintf(folder_name, "output/%02d_%02d_%02d_%02d_%02d_%02d", 1900 + ltm->tm_year,
          1 + ltm->tm_mon, ltm->tm_mday, ltm->tm_hour, ltm->tm_min,
          ltm->tm_sec);

  char command[100];
  sprintf(command, "mkdir output/%02d_%02d_%02d_%02d_%02d_%02d", 1900 + ltm->tm_year,
          1 + ltm->tm_mon, ltm->tm_mday, ltm->tm_hour, ltm->tm_min,
          ltm->tm_sec);

  system((const char *)command);

  if (system(folder_name) != 0) {
    debug_printf("Could not create output folder with name %s\n", folder_name);
    /*throw runtime_error("Could not create output folder with name " +
                        string(folder_name));*/
  }

  return string(folder_name);
}

void latest_outout(string folder){
  // check if output/latest folder exists in the current folder
  // if it does exist delete it
  if (system("ls | grep output/00_latest") > 0) {
    debug_print("Deleting latest output folder\n");
    system("rm -r output/00_latest");
  }

  // copy over the latest output folder to output/latest
  char command[100];
  sprintf(command, "cp -r %s output/00_latest", folder.c_str());
  system((const char *)command);
}

// TODO CHECK TYPES IF ALL CORRECT
// TODO arrays need to be initialized likely with 0 values
int main() {
  string folder_name = prepare_output(); // TODO delete empty folders

  debug_printf("Output folder: %s\n", folder_name.c_str());

  // Lattice Boltzmann Simulation in 2D
  debug_print("Starting\n");

  // TODO creating array directly segfaults because of size
  double ***F = malloc_3d(Ny, Nx, NL);
  debug_print("Initializing\n");

  srand(42);  // some seed

  // Initialize F
  for (int i = 0; i < Ny; i++) {
    for (int j = 0; j < Nx; j++) {
      for (int k = 0; k < NL; k++) {
        double rand_val = ((double) rand() / (RAND_MAX)) + 1;
        F[i][j][k] = 1 + 0.01 * rand_val;
      }
    }
  }

  int **x_coords = malloc_2d(Ny, Nx);
  int **y_coords = malloc_2d(Ny, Nx);

  meshgrid(x_coords, y_coords);

  // F[:,:,3] += 2 * (1+0.2*np.cos(2*np.pi*X/Nx*4))
  for (int i = 0; i < Ny; i++) {
    for (int j = 0; j < Nx; j++) {
      // in python we access [3] but here we do on [1] what this dictates is the direction we go
      // maybe something else too, but my brain is more fried than a kfc chicken
      F[i][j][1] += 2.0 * (1.0 + 0.2 * cos(2.0 * M_PI * (double)x_coords[i][j] / (double)Nx * 4.0));
    }
  }

  // rho = np.sum(F, axis=2)
  double **rho = malloc_2d_double(Ny, Nx);
  
  // reset rho
  for (int j = 0; j < Ny; j++) {
    for (int k = 0; k < Nx; k++) {
      rho[j][k] = 0;
    }
  }

  for (int i = 0; i < Ny; i++) {
    for (int j = 0; j < Nx; j++) {
      rho[i][j] = 0;
      for (int k = 0; k < NL; k++) {
        rho[i][j] += F[i][j][k];
      }
    }
  }

  // 	for i in idxs:		F[:,:,i] *= rho0 / rho
  for (int i = 0; i < NL; i++) {
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nx; k++) {
        F[j][k][i] *= rho0 / rho[j][k];
      }
    }
  }

  // Cylinder boundary
  int **cylinderX = malloc_2d(Ny, Nx);
  int **cylinderY = malloc_2d(Ny, Nx);
  double **cylinder = malloc_2d_double(Ny, Nx);

  // TODO we could reuse the meshgrid from above
  meshgrid(cylinderX, cylinderY);

  // cylinder = (X - Nx/4)**2 + (Y - Ny/2)**2 < (Ny/4)**2
  for (int i = 0; i < Ny; i++) {
    for (int j = 0; j < Nx; j++) {
      if ((pow((double)x_coords[i][j] - (double)Nx / 4, 2) + pow((double)y_coords[i][j] - (double)Ny / 2, 2)) <
          pow(Ny / 4, 2)) {
        cylinder[i][j] = 1;
      } else {
        cylinder[i][j] = 0;
      }
    }
  }
  

  // Simulation loop
  for (int i = 0; i < Nt; i++) {
    printf("Timestep %05d\n", i);

    //# Drift
    for (int j = 0; j < NL; j++) {

      double **temp = malloc_2d_double(Ny, Nx);

      for (int k = 0; k < Ny; k++) {
        for (int l = 0; l < Nx; l++) {
          temp[k][l] = F[k][l][j];
        }
      }

      roll2D(temp, Ny, Nx, cxs[j], 1);

      for (int k = 0; k < Ny; k++) {
        for (int l = 0; l < Nx; l++) {
          F[k][l][j] = temp[k][l];
        }
      }

      free(temp);

      temp = malloc_2d_double(Ny, Nx);
      
      for (int k = 0; k < Ny; k++) {
        for (int l = 0; l < Nx; l++) {
          temp[k][l] = F[k][l][j];
        }
      }

      roll2D(temp, Ny, Nx, cys[j], 0);

      for (int k = 0; k < Ny; k++) {
        for (int l = 0; l < Nx; l++) {
          F[k][l][j] = temp[k][l];
        }
      }

      free(temp);
      
    }


    // bndryF = F[cylinder,:]
    // its 2d of size 1941x9 but no idea how this is calculated ??????
    // TODO to support dynamic sized we could evaluate the size of the array and then allocate memory
    // for now its hardcoded :pikashrug:
    double **bndryF = malloc_2d_double(1941, NL);
    int index_bndryF = 0;
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nx; k++) {
        if (cylinder[j][k] == 1) {
          for (int l = 0; l < NL; l++) {
            bndryF[index_bndryF][l] = F[j][k][l];
          }
          index_bndryF++;
        }
      }
    }

                                        // 0,1,2,3,4,5,6,7,8  INDEXES
    // reorder columns bndryF = bndryF[:,[0,5,6,7,8,1,2,3,4]]
    for (int j = 0; j < 1941; j++) {
        // we love pointers dont we?
        double temp = bndryF[j][1];
        bndryF[j][1] = bndryF[j][5];
        bndryF[j][5] = temp;

        temp = bndryF[j][2];
        bndryF[j][2] = bndryF[j][6];
        bndryF[j][6] = temp;

        temp = bndryF[j][3];
        bndryF[j][3] = bndryF[j][7];
        bndryF[j][7] = temp;

        temp = bndryF[j][4];
        bndryF[j][4] = bndryF[j][8];
        bndryF[j][8] = temp;
    }

    // rho = np.sum(F,2)
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nx; k++) {
        rho[j][k] = 0;
        for (int l = 0; l < NL; l++) {
          rho[j][k] += F[j][k][l];
        }
      }
    }


    // ux = np.sum(F * cxs, 2) / rho
    double **ux = malloc_2d_double(Ny, Nx);
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nx; k++) {
        ux[j][k] = 0;
        for (int l = 0; l < NL; l++) {
          ux[j][k] += F[j][k][l] * cxs[l];
        }
        ux[j][k] /= rho[j][k];
      }
    }

    // uy = np.sum(F * cys, 2) / rho
    double **uy = malloc_2d_double(Ny, Nx);
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nx; k++) {
        uy[j][k] = 0;
        for (int l = 0; l < NL; l++) {
          uy[j][k] += F[j][k][l] * cys[l];
        }
        uy[j][k] /= rho[j][k];
      }
    }

    double ***Feq = malloc_3d(Ny, Nx, NL);
    // set to zero
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nx; k++) {
        for (int l = 0; l < NL; l++) {
          Feq[j][k][l] = 0;
        }
      }
    }

    for (int k = 0; k < NL; k++) {
      for (int j = 0; j < Ny; j++) {
        for (int l = 0; l < Nx; l++) {

          double rho_val = rho[j][l];
          double weight_val = weights[k];

          // 3*(cx*ux+cy*uy)
          double first = 3 * (cxs[k] * ux[j][l] + cys[k] * uy[j][l]);

          // 9*(cx*ux+cy*uy)**2/2
          double second = 9 * pow(cxs[k] * ux[j][l] + cys[k] * uy[j][l], 2) / 2;

          // 3*(ux**2+uy**2)/2
          double third = 3 * (pow(ux[j][l], 2) + pow(uy[j][l], 2)) / 2;

          Feq[j][l][k] = rho_val * weight_val * (1 + first + second - third);
        }
      }
    }

    // F += -(1.0/tau) * (F - Feq)
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nx; k++) {
        for (int l = 0; l < NL; l++) {
          F[j][k][l] += -(1.0 / tau) * (F[j][k][l] - Feq[j][k][l]);
        }
      }
    }

    debug_print("F calculated\n");

    // Apply boundary
    // F[cylinder,:] = bndryF
    int index_bndryF2 = 0;
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nx; k++) {
        if (cylinder[j][k] == 1) {
          for (int l = 0; l < NL; l++) {
            F[j][k][l] = bndryF[index_bndryF2][l];
          }
          index_bndryF2++;
        }
      }
    }

    // set ux and uy to zero where cylinder is 1
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nx; k++) {
        if (cylinder[j][k] == 1) {
          ux[j][k] = 0;
          uy[j][k] = 0;
        }
      }
    }

    // vorticity = (np.roll(ux, -1, axis=0) - np.roll(ux, 1, axis=0)) -
    // (np.roll(uy, -1, axis=1) - np.roll(uy, 1, axis=1)) vorticity[cylinder]
    // = np.nan vorticity = np.ma.array(vorticity, mask=cylinder)
    double **vorticity = malloc_2d_double(Ny, Nx);
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nx; k++) {

        // (np.roll(ux, -1, axis=0) - np.roll(ux, 1, axis=0))
        double ux_roll = ux[j][(k - 1 + Nx) % Nx] - ux[j][(k + 1 + Nx) % Nx];

        // (np.roll(uy, -1, axis=1) - np.roll(uy, 1, axis=1))
        double uy_roll = uy[(j - 1 + Ny) % Ny][k] - uy[(j + 1 + Ny) % Ny][k];

        vorticity[j][k] = ux_roll - uy_roll;

        if (cylinder[j][k] == 1) {
          vorticity[j][k] = 0; 
        }
      }
    }

    char vortex_filename[100];
    sprintf(vortex_filename, "%s/vorticity_%05d.npy", folder_name.c_str(), i);
    // TODO for benchmarking only save vorticity from the last step
    save_npy_2d_double(vorticity, Ny, Nx, vortex_filename);
  }

  latest_outout(folder_name);

  return 0;
}

// to run and compile run: clear && g++ -o 2D/2D 2D/main.cpp -lm && 2D/2D
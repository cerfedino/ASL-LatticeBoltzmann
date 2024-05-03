#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <cmath>
#include <exception>
#include <stdexcept>
#include <vector>
#include <string>
   
#include "include/npy.hpp"
#include "include/utils.h"

#ifdef DEBUG
#define debug_printf(fmt, ...)  fprintf(stdout, fmt, __VA_ARGS__)
#define debug_print(fmt)        fprintf(stdout, fmt)

#else
// If DEBUG is not defined we expand macros to whitespace
#define debug_printf(fmt, ...);
#define debug_print(fmt);
#define save_npy_3d_double(array, x, y, z, filename);
#define save_npy_2d_double(array, x, y, filename); 
#define save_npy_2d_int(array, x, y, filename); 
#define make_output_folder() ""
#define make_latest_output(folder);
#endif

using namespace std;

int Nx = 400;   // resolution in x
int Ny = 100;   // resolution in y
#define rho0 0.01 // reciprocal average density
#define tau -1.66666667  // reciprocal collision timescale (1/0.6)
#define Nt 5000  // number of timesteps

// #define Nt 30  // number of timesteps

// Lattice speeds / weights
#define NL 9
const double idx[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
const double cxs[9] = {0.0, 0.0, 1.0, 1.0, 1.0, 0.0, -1.0, -1.0, -1.0};
const double cys[9] = {0.0, 1.0, 1.0, 0.0, -1.0, -1.0, -1.0, 0.0, 1.0};
const double weights[9] = {4.0 / 9, 1.0 / 9,  1.0 / 36, 1.0 / 9, 1.0 / 36,
                     1.0 / 9, 1.0 / 36, 1.0 / 9,  1.0 / 36};  // sums to 1


void meshgrid(int *x_coords, int *y_coords) {
  for (int i = 0; i < Ny; ++i) {
    for (int j = 0; j < Nx; ++j) {
      x_coords[i*Nx +j] = j;
      y_coords[i*Nx +j] = i;
    }
  }
}


// TODO CHECK TYPES IF ALL CORRECT
// TODO arrays need to be initialized likely with 0 values
inline int run() {
  string folder_name = make_output_folder(); // TODO delete empty folders

  debug_printf("Output folder: %s\n", folder_name.c_str());

  // Lattice Boltzmann Simulation in 2D
  debug_print("Starting\n");

  // double *BIG_CHUNGUS = (double *) malloc((2*(Ny*Nx*NL)+6*(Ny*Nx)+(1941*NL)) * sizeof(double) + 2*(Ny*Nx)*sizeof(int));
  // debug_printf("Allocating %lu bytes for BIG_CHUNGUS\n", (2*(Ny*Nx*NL)+5*(Ny*Nx)+(1941*NL)) * sizeof(double) + 2*(Ny*Nx)*sizeof(int));
  // if (BIG_CHUNGUS == NULL) {
  //   printf("Memory allocation failed for BIG_CHUNGUS\n");
  //   return 1;
  // }

  

  double *Feq =       (double *) malloc(Ny*Nx*NL*sizeof(double));
  double *F =         (double *) malloc(Ny*Nx*NL*sizeof(double));
  double *vorticity = (double *) malloc(Ny*Nx*sizeof(double));
  double *rho =       (double *) calloc(Ny*Nx,sizeof(double));
  double *cylinder =  (double *) malloc(Ny*Nx*sizeof(double));
  double *ux =        (double *) malloc(Ny*Nx*sizeof(double));
  double *uy =        (double *) malloc(Ny*Nx*sizeof(double)); 
  double *temp =      (double *) malloc(Ny*Nx*sizeof(double));
  double *bndryF =    (double *) malloc(1941*NL*sizeof(double));
  int *x_coords =        (int *) malloc(Ny*Nx*sizeof(int));
  int *y_coords =        (int *) malloc(Ny*Nx*sizeof(int));


  debug_print("Initializing\n");

  srand(42);  // some seed

  // Initialize F
  for (int i = 0; i < Ny; i++) {
    for (int j = 0; j < Nx; j++) {
      for (int k = 0; k < NL; k++) {
        double rand_val = ((double) rand() / (RAND_MAX)) + 1;
        F[i*(Nx*NL)+ j*NL + k] = 1 + 0.01 * rand_val;
      }
    }
  }



  meshgrid(x_coords, y_coords);

  // F[:,:,3] += 2 * (1+0.2*np.cos(2*np.pi*X/Nx*4))
  for (int i = 0; i < Ny; i++) {
    for (int j = 0; j < Nx; j++) {
      // in python we access [3] but here we do on [1] what this dictates is the direction we go
      // maybe something else too, but my brain is more fried than a kfc chicken
      // AHAHAHAHAHAHAH HILARIOUS KARLO LOL IM LITERALLY DYING OF LAUGHTER
      F[i*(Nx*NL) +j*NL +1] += 2.0 * (1.0 + 0.2 * cos(2.0 * M_PI * (double)x_coords[i*Nx +j] / (double)Nx * 4.0));
    }
  }


  int res = 0;
  for (int i = 0; i < Ny; i++) {
    for (int j = 0; j < Nx; j++) {
      res = 0;
      for (int k = 0; k < NL; k++) {
        res += F[i*(Nx*NL)+ j*NL +k];
      }
      rho[i*Nx +j] = res;
    }
  }

  // 	for i in idxs:		F[:,:,i] *= rho0 / rho
  for (int j = 0; j < Ny; j++) {
    for (int k = 0; k < Nx; k++) {
      for (int i = 0; i < NL; i++) {
        F[j*(Nx*NL)+ k*NL +i] *= rho0 * rho[j*Nx +k];
      }
    }
  }



  // cylinder = (X - Nx/4)**2 + (Y - Ny/2)**2 < (Ny/4)**2
  for (int i = 0; i < Ny; i++) {
    for (int j = 0; j < Nx; j++) {
      cylinder[i*Nx +j] = (pow((double)x_coords[i*Nx +j] - (double)Nx / 4, 2) + pow((double)y_coords[i*Nx +j] - (double)Ny / 2, 2)) < pow(Ny / 4, 2);
    }
  }
  
  
  // Simulation loop
  for (int i = 0; i < Nt; i++) {
    debug_printf("\r%d", i);

    //# Drift
    for (int j = 0; j < NL; j++) {

      // Roll 
      int shiftX = cys[j];
      int shiftY = cxs[j];
      for (int k = 0; k < Ny; k++) {
        for (int l = 0; l < Nx; l++) {
          temp[((k + shiftY + Ny) % Ny)*Nx +((l + shiftX + Nx) % Nx)] = F[k*(Nx*NL)+ l*NL +j];
        }
      }

      for (int k = 0; k < Ny; k++) {
        for (int l = 0; l < Nx; l++) {
          F[k*(Nx*NL)+ l*NL +j] = temp[k*Nx +l];
        }
      }
    }


    // bndryF = F[cylinder,:]
    // its 2d of size 1941x9 but no idea how this is calculated ??????
    // TODO to support dynamic sized we could evaluate the size of the array and then allocate memory
    // for now its hardcoded :pikashrug:
    // TODO: FIX
    int index_bndryF = 0;
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nx; k++) {
        if (cylinder[j*Nx +k] == 1) {
          for (int l = 0; l < NL; l++) {
            bndryF[index_bndryF*NL +l] = F[j*(Nx*NL)+ k*NL +l];
          }
          index_bndryF++;
        }
      }
    }

                                        // 0,1,2,3,4,5,6,7,8  INDEXES
    // reorder columns bndryF = bndryF[:,[0,5,6,7,8,1,2,3,4]]
    for (int j = 0; j < 1941; j++) {
        // we love pointers dont we?
        double temp = bndryF[j*NL +1];
        bndryF[j*NL +1] = bndryF[j*NL +5];
        bndryF[j*NL +5] = temp;

        temp = bndryF[j*NL +2];
        bndryF[j*NL +2] = bndryF[j*NL +6];
        bndryF[j*NL +6] = temp;

        temp = bndryF[j*NL +3];
        bndryF[j*NL +3] = bndryF[j*NL +7];
        bndryF[j*NL +7] = temp;

        temp = bndryF[j*NL +4];
        bndryF[j*NL +4] = bndryF[j*NL +8];
        bndryF[j*NL +8] = temp;
    }

    // rho = np.sum(F,2)
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nx; k++) {
        double res1 = 0;
        double res2 = 0;
        double res3 = 0;
        for (int l = 0; l < NL; l++) {
          res1 += F[j*(Nx*NL)+ k*NL +l];
          res2 += F[j*(Nx*NL)+ k*NL +l] * cxs[l];
          res3 += F[j*(Nx*NL)+ k*NL +l] * cys[l];
        }
        rho[j*Nx +k] = res1;
        ux[j*Nx +k] = res2 / res1;
        uy[j*Nx +k] = res3 / res1;
      }
    }

    // set to zero
    memset(Feq, 0, Ny*Nx*NL*sizeof(double));

    for (int k = 0; k < NL; k++) {
      for (int j = 0; j < Ny; j++) {
        for (int l = 0; l < Nx; l++) {

          double rho_val = rho[j*Nx +l];
          double weight_val = weights[k];

          // 3*(cx*ux+cy*uy)
          double first = 3 * (cxs[k] * ux[j*Nx +l] + cys[k] * uy[j*Nx +l]);

          // 9*(cx*ux+cy*uy)**2/2
          double second = 9 * pow(cxs[k] * ux[j*Nx +l] + cys[k] * uy[j*Nx +l], 2) / 2;

          // 3*(ux**2+uy**2)/2
          double third = 3 * (pow(ux[j*Nx +l], 2) + pow(uy[j*Nx +l], 2)) / 2;

          Feq[j*(Nx*NL)+ l*NL +k] = rho_val * weight_val * (1 + first + second - third);
        }
      }
    }

    // F += -(1.0/tau) * (F - Feq)
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nx; k++) {
        for (int l = 0; l < NL; l++) {
          F[j*(Nx*NL)+ k*NL +l] += tau * (F[j*(Nx*NL)+ k*NL +l] - Feq[j*(Nx*NL)+ k*NL +l]);
        }
      }
    }


    // Apply boundary
    // F[cylinder,:] = bndryF
    int index_bndryF2 = 0;
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nx; k++) {
        if (cylinder[j*Nx +k] == 1) {
          for (int l = 0; l < NL; l++) {
            F[j*(Nx*NL)+ k*NL +l] = bndryF[index_bndryF2*NL +l];
          }
          index_bndryF2++;
        }
      }
    }

    

    // set ux and uy to zero where cylinder is 1
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nx; k++) {
        if (cylinder[j*Nx +k] == 1) {
          ux[j*Nx +k] = 0;
          uy[j*Nx +k] = 0;
        }
      }
    }

    // vorticity = (np.roll(ux, -1, axis=0) - np.roll(ux, 1, axis=0)) -
    // (np.roll(uy, -1, axis=1) - np.roll(uy, 1, axis=1)) vorticity[cylinder]
    // = np.nan vorticity = np.ma.array(vorticity, mask=cylinder)
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nx; k++) {

        // (np.roll(ux, -1, axis=0) - np.roll(ux, 1, axis=0))
        double ux_roll = ux[j*Nx +(k - 1 + Nx) % Nx] - ux[j*Nx +(k + 1 + Nx) % Nx];

        // (np.roll(uy, -1, axis=1) - np.roll(uy, 1, axis=1))
        double uy_roll = uy[((j - 1 + Ny) % Ny)*Nx +k] - uy[((j + 1 + Ny) % Ny)*Nx +k];

        vorticity[j*Nx +k] = cylinder[j*Nx +k] == 1 ? 0 : ux_roll - uy_roll;
      }
    }

    #ifdef DEBUG
    char vortex_filename[100];
    sprintf(vortex_filename, "%s/vorticity_%05d.npy", folder_name.c_str(), i);
    // TODO for benchmarking only save vorticity from the last step
    save_npy_2d_double(vorticity, Ny, Nx, vortex_filename);
    #endif

    
    
    
    
    
  }
  printf("\n");

  free(ux);
  free(uy);
  free(Feq);
  free(vorticity);
  free(bndryF);
  free(F);
  free(x_coords);
  free(y_coords);
  free(rho);
  free(cylinder);
  free(temp);
  // free(BIG_CHUNGUS);
  make_latest_output(folder_name);

  return 0;
}




int main(int argc, char const *argv[])
{
  std::cin>>Nx>>Ny;

  unsigned long long start_cycle, end_cycle;
  time_t start_sec, end_sec;
  asm volatile ("RDTSC" : "=A" (start_cycle));
  time(&start_sec);
  run();
  asm volatile ("RDTSC" : "=A" (end_cycle));
  time(&end_sec);
  printf("Cycles taken: %llu (%ld seconds)\n", end_cycle - start_cycle, end_sec - start_sec);

  return 0;
}

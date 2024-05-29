// description = "....."
// flops = Ny*Nx*NL(1 div + 1 add + 1 add + 1 mult)+ Ny*Nx( 2 add + 5 mult + 1
// cos + 1 div )+Ny*Nx*NL*(1 add)+Ny*Nx*NL*(2 adds)+ Ny*Nx*(2 pow + 2 subs + 2
// divs+ 1 add)+ Nt*(Ny*Nx*NL(3 adds + 2 mults) + Ny*Nx*(2 divs))+NL*Ny*Nx*(9
// mults + 2 div  + 3 pows + 6 adds)+Ny*Nx*NL (2 add + 1 mult )+Ny*Nx*(3 adds)
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <cmath>
#include <exception>
#include <immintrin.h>
#include <stdexcept>
#include <string>
#include <vector>

#include "include/npy.hpp"
#include "include/profile.h"
#include "include/utils.c"

#ifdef DEBUG
#define debug_printf(fmt, ...) fprintf(stdout, fmt, __VA_ARGS__)
#define debug_print(fmt) fprintf(stdout, fmt)

#else
// If DEBUG is not defined we expand macros to whitespace
#define debug_printf(fmt, ...) ;
#define debug_print(fmt) ;
#define save_npy_3d_double(array, x, y, z, filename) ;
#define save_npy_2d_double(array, x, y, filename) ;
#define save_npy_2d_int(array, x, y, filename) ;
#define make_output_folder() ""
#define make_latest_output(folder) ;
#endif

using namespace std;



#ifdef MNx
#define Nx (MNx)
#endif
#ifdef MNy
#define Ny (MNy)
#endif
#ifdef MNt
#define Nt (MNt)
#endif


#define rho0 0.01       // reciprocal average density
#define tau -1.66666667 // reciprocal collision timescale (1/0.6)
#define tau_plus_1 tau + 1

// #define Nt 30  // number of timesteps

// Lattice speeds / weights
#define NL 9
const double idx[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
const double cxs[9] = {0.0, 0.0, 1.0, 1.0, 1.0, 0.0, -1.0, -1.0, -1.0};
const double cys[9] = {0.0, 1.0, 1.0, 0.0, -1.0, -1.0, -1.0, 0.0, 1.0};
const double weights[9] = {tau * 4.0 / 9, tau * 1.0 / 9, tau * 1.0 / 36, tau * 1.0 / 9, tau * 1.0 / 36, tau * 1.0 / 9, tau * 1.0 / 36, tau * 1.0 / 9, tau * 1.0 / 36}; // sums to tau * 1

void meshgrid(int *x_coords, int *y_coords) {
  for (int i = 0; i < Ny; ++i) {
    for (int j = 0; j < Nx; ++j) {
      x_coords[i * Nx + j] = j;
      y_coords[i * Nx + j] = i;
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

  // double *BIG_CHUNGUS = (double *) malloc((2*(Ny*Nx*NL)+6*(Ny*Nx)+(1941*NL))
  // * sizeof(double) + 2*(Ny*Nx)*sizeof(int)); debug_printf("Allocating %lu
  // bytes for BIG_CHUNGUS\n", (2*(Ny*Nx*NL)+5*(Ny*Nx)+(1941*NL)) *
  // sizeof(double) + 2*(Ny*Nx)*sizeof(int)); if (BIG_CHUNGUS == NULL) {
  //   printf("Memory allocation failed for BIG_CHUNGUS\n");
  //   return 1;
  // }

  double *Feq = (double *)aligned_alloc(32, Ny * Nx * NL * sizeof(double));
  double *F = (double *)aligned_alloc(32, Ny * Nx * NL * sizeof(double));
  double *vorticity = (double *)aligned_alloc(32, Ny * Nx * sizeof(double));
  double *rho = (double *)aligned_alloc(32, Ny * Nx * sizeof(double));
  double *cylinder = (double *)aligned_alloc(32, Ny * Nx * sizeof(double));
  double *ux = (double *)aligned_alloc(32, Ny * Nx * sizeof(double));
  double *uy = (double *)aligned_alloc(32, Ny * Nx * sizeof(double));
  double *temp = (double *)aligned_alloc(32, Ny * Nx * sizeof(double));
  int *x_coords = (int *)aligned_alloc(32, Ny * Nx * sizeof(int));
  int *y_coords = (int *)aligned_alloc(32, Ny * Nx * sizeof(int));

  debug_print("Initializing\n");

  srand(42); // some seed

  // Initialize F
  // flops = Ny*Nx*NL(1 div + 1 add + 1 add + 1 mult)
  for (int i = 0; i < Ny; i++) {
    for (int j = 0; j < Nx; j++) {
      for (int k = 0; k < NL; k++) {
        double rand_val = ((double)rand() / (RAND_MAX)) + 1;
        F[i * (Nx * NL) + k * Nx + j] = 1 + 0.01 * rand_val;
      }
    }
  }

  meshgrid(x_coords, y_coords);

  // F[:,:,3] += 2 * (1+0.2*np.cos(2*np.pi*X/Nx*4))
  // flops = Ny*Nx( 2 add + 5 mult + 1 cos + 1 div )
  for (int i = 0; i < Ny; i++) {
    for (int j = 0; j < Nx; j++) {
      // in python we access [3] but here we do on [1] what this dictates is the
      // direction we go maybe something else too, but my brain is more fried
      // than a kfc chicken AHAHAHAHAHAHAH HILARIOUS KARLO LOL IM LITERALLY
      // DYING OF LAUGHTER
      F[i * (Nx * NL) + Nx + j] += 2.0 * (1.0 + 0.2 * cos(2.0 * M_PI * (double)x_coords[i * Nx + j] / (double)Nx * 4.0));
    }
  }

  // flops = Ny*Nx*NL*(1 add)
  int res = 0;
  for (int i = 0; i < Ny; i++) {
    for (int j = 0; j < Nx; j++) {
      res = 0;
      for (int k = 0; k < NL; k++) {
        res += F[i * (Nx * NL) + k * Nx + j];
      }
      rho[i * Nx + j] = res;
    }
  }

  // 	for i in idxs:		F[:,:,i] *= rho0 / rho
  // flops = Ny*Nx*NL*(2 adds)
  for (int j = 0; j < Ny; j++) {
    for (int k = 0; k < Nx; k++) {
      for (int i = 0; i < NL; i++) {
        F[j * (Nx * NL) + i * Nx + k] *= rho0 * rho[j * Nx + k];
      }
    }
  }

  // cylinder = (X - Nx/4)**2 + (Y - Ny/2)**2 < (Ny/4)**2
  // flops = Ny*Nx*(2 pow + 2 subs + 2 divs+ 1 add)
  for (int i = 0; i < Ny; i++) {
    for (int j = 0; j < Nx; j++) {
      cylinder[i * Nx + j] = (pow((double)x_coords[i * Nx + j] - (double)Nx / 4, 2) + pow((double)y_coords[i * Nx + j] - (double)Ny / 2, 2)) < pow(Ny / 4, 2);
    }
  }

  int bndryF_size = 0;

  // loop trough cylinder and count the number of 1s
  // flops = 0 just reding the values
  for (int i = 0; i < Ny; i++) {
    for (int j = 0; j < Nx; j++) {
      if (cylinder[i * Nx + j] == 1) {
        bndryF_size++;
      }
    }
  }

  double *bndryF = (double *)aligned_alloc(32, bndryF_size * NL * sizeof(double));

    // Currently counting only compulsory misses
  profiler *rho_profiler = init_profiler(5 * Ny * Nx * NL + 2 * Ny * Nx, 8 * 5 * Ny * Nx * NL + 3 * Ny * Nx), *feq_profiler = init_profiler(19 * Nx * Ny * NL, 8 * 13 * Nx * Ny * NL), *f_profiler = init_profiler(2 * Nx * Ny * NL, 8 * 3 * Nx * Ny * NL),
           *vort_profiler = init_profiler(2 * Nx * Ny, 8 * 6 * Nx * Ny);

  // flops = Nt*(Ny*Nx*NL(3 adds + 2 mults) + Ny*Nx*(2 divs))+NL*Ny*Nx*(9 mults
  // + 2 div  + 3 pows + 6 adds)+Ny*Nx*NL (2 add + 1 mult )+Ny*Nx*(3 adds)
#ifdef PROFILE
  const int PROFILE_RUNS = 5;
  const int PROFILE_DIGITS = floor(log10(PROFILE_RUNS)) + 1;
  printf("\rRun %-*d/%d done", PROFILE_DIGITS, 0, PROFILE_RUNS);
  fflush(stdout);
  for (int profile_counter = 0; profile_counter < PROFILE_RUNS; profile_counter++) {
#endif
    // Simulation loop
    for (int i = 0; i < Nt; i++) {
      debug_printf("\r%d", i);

      // # Drift
      // flops = 0
      for (int j = 0; j < NL; j++) {

        // Roll
        int shiftX = cys[j];
        int shiftY = cxs[j];
        for (int k = 0; k < Ny; k++) {
          for (int l = 0; l < Nx; l++) {
            temp[((k + shiftY + Ny) % Ny) * Nx + ((l + shiftX + Nx) % Nx)] = F[k * (Nx * NL) + j * Nx + l];
          }
        }

        for (int k = 0; k < Ny; k++) {
          for (int l = 0; l < Nx; l++) {
            F[k * (Nx * NL) + j * Nx + l] = temp[k * Nx + l];
          }
        }
      }

      // bndryF = F[cylinder,:]
      // its 2d of size 1941x9 but no idea how this is calculated ??????
      // TODO to support dynamic sized we could evaluate the size of the
      // array and then allocate memory for now its hardcoded :pikashrug:
      // TODO: FIX
      // flops = 0
      int index_bndryF = 0;
      for (int j = 0; j < Ny; j++) {
        for (int k = 0; k < Nx; k++) {
          if (cylinder[j * Nx + k] == 1) {
            for (int l = 0; l < NL; l++) {
              bndryF[index_bndryF * NL + l] = F[j * (Nx * NL) + l * Nx + k];
            }
            index_bndryF++;
          }
        }
      }

      // 0,1,2,3,4,5,6,7,8  INDEXES
      // reorder columns bndryF = bndryF[:,[0,5,6,7,8,1,2,3,4]]
      // flops = 0
      for (int j = 0; j < bndryF_size; j++) {
        // we love pointers dont we?
        double temp = bndryF[j * NL + 1];
        bndryF[j * NL + 1] = bndryF[j * NL + 5];
        bndryF[j * NL + 5] = temp;

        temp = bndryF[j * NL + 2];
        bndryF[j * NL + 2] = bndryF[j * NL + 6];
        bndryF[j * NL + 6] = temp;

        temp = bndryF[j * NL + 3];
        bndryF[j * NL + 3] = bndryF[j * NL + 7];
        bndryF[j * NL + 7] = temp;

        temp = bndryF[j * NL + 4];
        bndryF[j * NL + 4] = bndryF[j * NL + 8];
        bndryF[j * NL + 8] = temp;
      }

      // rho = np.sum(F,2)
      // flops = Ny*Nx*NL(3 adds + 2 mults) + Ny*Nx*(2 divs)
      start_run(rho_profiler);
      for (int j = 0; j < Ny; j++) {
        for (int k = 0; k < Nx; k++) {
          double res1 = 0;
          double res2 = 0;
          double res3 = 0;
          for (int l = 0; l < NL; l++) {
            res1 += F[j * (Nx * NL) + l * Nx + k];
            res2 += F[j * (Nx * NL) + l * Nx + k] * cxs[l];
            res3 += F[j * (Nx * NL) + l * Nx + k] * cys[l];
          }
          double inv = 1 / res1;
          rho[j * Nx + k] = res1;
          ux[j * Nx + k] = res2 * inv;
          uy[j * Nx + k] = res3 * inv;
        }
      }
      end_run(rho_profiler);

      // set to zero
      memset(Feq, 0, Ny * Nx * NL * sizeof(double));

      // flops = NL*Ny*Nx*(9 mults + 2 div  + 3 pows + 6 adds)
      __m256d const_vec_1 = _mm256_set1_pd(1);
      __m256d const_vec_3 = _mm256_set1_pd(3);
      __m256d const_vec_4_5 = _mm256_set1_pd(4.5);
      __m256d const_vec_1_5 = _mm256_set1_pd(1.5);
      start_run(feq_profiler);
      for (int k = 0; k < NL; k++) {
        __m256d x_vec = _mm256_set1_pd(cxs[k]);
        __m256d y_vec = _mm256_set1_pd(cys[k]);
        __m256d weights_vec = _mm256_set1_pd(weights[k]);
        for (int j = 0; j < Ny; j++) {
          int l;
          for (l = 0; l < Nx - 3; l += 4) {
            __m256d rho_vec = _mm256_load_pd(&rho[j * Nx + l]); // rho[j * Nx + l]

            __m256d ux_vec = _mm256_load_pd(&ux[j * Nx + l]);             // ux[j * Nx + l]
            __m256d uy_vec = _mm256_load_pd(&uy[j * Nx + l]);             // uy[j * Nx + l]
            __m256d ux_mul_vec = _mm256_mul_pd(x_vec, ux_vec);            // cxs[k] * ux[j * Nx + l]
            __m256d uy_mul_vec = _mm256_mul_pd(y_vec, uy_vec);            // cys[k] * uy[j * Nx + l]
            __m256d uxuy_add_vec = _mm256_add_pd(ux_mul_vec, uy_mul_vec); // cxs[k] * ux[j * Nx + l] + cys[k] * uy[j * Nx + l]

            __m256d first_vec = _mm256_mul_pd(uxuy_add_vec, const_vec_3); // 3 * (cxs[k] * ux[j * Nx + l] + cys[k] * uy[j * Nx + l])

            __m256d uxuypow_vec = _mm256_mul_pd(uxuy_add_vec, uxuy_add_vec); // pow(cxs[k] * ux[j * Nx + l] + cys[k] * uy[j * Nx + l], 2)
            __m256d second_vec = _mm256_mul_pd(uxuypow_vec, const_vec_4_5);  // 9 * pow(cxs[k] * ux[j * Nx + l] + cys[k] * uy[j * Nx + l], 2) / 2

            __m256d uxpow_vec = _mm256_mul_pd(ux_vec, ux_vec);                // pow(ux[j * Nx + l], 2)
            __m256d uypow_vec = _mm256_mul_pd(uy_vec, uy_vec);                // pow(uy[j * Nx + l], 2)
            __m256d uxuypowadd_vec = _mm256_add_pd(uxpow_vec, uypow_vec);     // pow(ux[j * Nx + l], 2) + pow(uy[j * Nx + l], 2)
            __m256d third_vec = _mm256_mul_pd(uxuypowadd_vec, const_vec_1_5); // 3 * (pow(ux[j * Nx + l], 2) + pow(uy[j * Nx + l], 2)) / 2

            __m256d first_sum = _mm256_add_pd(const_vec_1, first_vec); // 1 + first
            __m256d second_sum = _mm256_sub_pd(second_vec, third_vec); // second - third
            __m256d first_mul = _mm256_mul_pd(rho_vec, weights_vec);   // rho_val * weight_val
            __m256d third_sum = _mm256_add_pd(first_sum, second_sum);  // 1 + first + second - third
            __m256d second_mul = _mm256_mul_pd(first_mul, third_sum);  // rho_val * weight_val * (1 + first + second - third)

            _mm256_store_pd(&Feq[j * (Nx * NL) + k * Nx + l], second_mul);
          }
          for (; l < Nx; l++) {
            double rho_val = rho[j * Nx + l];
            double weight_val = weights[k];

            // 3*(cx*ux+cy*uy)
            double first = 3 * (cxs[k] * ux[j * Nx + l] + cys[k] * uy[j * Nx + l]);

            // 9*(cx*ux+cy*uy)**2/2
            double second = 9 * pow(cxs[k] * ux[j * Nx + l] + cys[k] * uy[j * Nx + l], 2) / 2;

            // 3*(ux**2+uy**2)/2
            double third = 3 * (pow(ux[j * Nx + l], 2) + pow(uy[j * Nx + l], 2)) / 2;

            Feq[j * (Nx * NL) + k * Nx + l] = rho_val * weight_val * (1 + first + second - third);
          }
        }
      }
      end_run(feq_profiler);

      // F += -(1.0/tau) * (F - Feq)
      // flops = Ny*Nx*NL (1 add + 1 mult )
      __m256d const_vec_tau = _mm256_set1_pd(tau_plus_1);
      start_run(f_profiler);
      for (int j = 0; j < Ny; j++) {
        for (int l = 0; l < NL; l++) {
          int k;
          for (k = 0; k < Nx; k += 8) {
            __m256d feq_vec = _mm256_load_pd(Feq + j * (Nx * NL) + l * Nx + k);
            __m256d f_vec = _mm256_load_pd(F + j * (Nx * NL) + l * Nx + k);
            __m256d feq1_vec = _mm256_load_pd(Feq + j * (Nx * NL) + l * Nx + k + 4);
            __m256d f1_vec = _mm256_load_pd(F + j * (Nx * NL) + l * Nx + k + 4);

            __m256d res_vec = _mm256_fmsub_pd(const_vec_tau, f_vec, feq_vec);
            _mm256_store_pd(F + j * (Nx * NL) + l * Nx + k, res_vec);

            __m256d res_vec1 = _mm256_fmsub_pd(const_vec_tau, f1_vec, feq1_vec);
            _mm256_store_pd(F + j * (Nx * NL) + l * Nx + k + 4, res_vec1);
          }
          for (; k < Nx; k++) {
            F[j * (Nx * NL) + l * Nx + k] = tau_plus_1 * F[j * (Nx * NL) + l * Nx + k] - Feq[j * (Nx * NL) + l * Nx + k];
          }
        }
      }
      end_run(f_profiler);

      // Apply boundary
      // F[cylinder,:] = bndryF
      // flops = 0
      int index_bndryF2 = 0;
      for (int j = 0; j < Ny; j++) {
        for (int k = 0; k < Nx; k++) {
          if (cylinder[j * Nx + k] == 1) {
            for (int l = 0; l < NL; l++) {
              F[j * (Nx * NL) + l * Nx + k] = bndryF[index_bndryF2 * NL + l];
            }
            index_bndryF2++;
          }
        }
      }

      // set ux and uy to zero where cylinder is 1
      // flops = 0
      for (int j = 0; j < Ny; j++) {
        for (int k = 0; k < Nx; k++) {
          if (cylinder[j * Nx + k] == 1) {
            ux[j * Nx + k] = 0;
            uy[j * Nx + k] = 0;
          }
        }
      }

      // vorticity = (np.roll(ux, -1, axis=0) - np.roll(ux, 1, axis=0)) -
      // (np.roll(uy, -1, axis=1) - np.roll(uy, 1, axis=1)) vorticity[cylinder]
      // = np.nan vorticity = np.ma.array(vorticity, mask=cylinder)

      // flops = Ny*Nx*(3 adds)
      start_run(vort_profiler);
      for (int j = 0; j < Ny; j++) {
        // Calculate k = 0 boundary
        double ux_roll = ux[j * Nx + Nx - 1] - ux[j * Nx + 1];
        double uy_roll = uy[((j - 1 + Ny) % Ny) * Nx] - uy[((j + 1 + Ny) % Ny) * Nx];
        vorticity[j * Nx] = cylinder[j * Nx] == 1 ? 0 : ux_roll - uy_roll;

        // Calculate k = [1, Nx - 2]
        for (int k = 1; k < Nx - 1; k++) {
          // (np.roll(ux, -1, axis=0) - np.roll(ux, 1, axis=0))
          double ux_roll = ux[j * Nx + k - 1] - ux[j * Nx + k + 1];

          // (np.roll(uy, -1, axis=1) - np.roll(uy, 1, axis=1))
          double uy_roll = uy[((j - 1 + Ny) % Ny) * Nx + k] - uy[((j + 1) % Ny) * Nx + k];

          vorticity[j * Nx + k] = cylinder[j * Nx + k] == 1 ? 0 : ux_roll - uy_roll;
        }

        // Calculate k = Nx - 1 boundary
        ux_roll = ux[j * Nx + Nx - 2] - ux[j * Nx];
        uy_roll = uy[((j - 1 + Ny) % Ny) * Nx + Nx - 1] - uy[((j + 1) % Ny) * Nx + Nx - 1];
        vorticity[j * Nx + Nx - 1] = cylinder[j * Nx + Nx - 1] == 1 ? 0 : ux_roll - uy_roll;
      }
      end_run(vort_profiler);

#ifdef DEBUG
      char vortex_filename[100];
      sprintf(vortex_filename, "%s/vorticity_%05d.npy", folder_name.c_str(), i);
      // TODO for benchmarking only save vorticity from the last step
      save_npy_2d_double(vorticity, Ny, Nx, vortex_filename);
#endif
    }
#ifdef PROFILE
    printf("\rRun %-*d/%d done", PROFILE_DIGITS, profile_counter + 1, PROFILE_RUNS);
    fflush(stdout);
  }
#endif
#ifdef PROFILE
  printf("\nProfiling results:\n");
  profiler_stats rho_stats = finish_profiler(rho_profiler);
  printf("- Rho  Calculation: %4.2f Flops/Cycle, %10ld cycles in %d runs. "
         "Arithmetic intensity: %4.2f\n",
         rho_stats.performance, rho_stats.cycles, rho_stats.runs, rho_stats.arithmetic_intensity);
  profiler_stats feq_stats = finish_profiler(feq_profiler);
  printf("- FEQ  Calculation: %4.2f Flops/Cycle, %10ld cycles in %d runs. "
         "Arithmetic intensity: %4.2f\n",
         feq_stats.performance, feq_stats.cycles, feq_stats.runs, feq_stats.arithmetic_intensity);
  profiler_stats f_stats = finish_profiler(f_profiler);
  printf("- F    Calculation: %4.2f Flops/Cycle, %10ld cycles in %d runs. "
         "Arithmetic intensity: %4.2f\n",
         f_stats.performance, f_stats.cycles, f_stats.runs, f_stats.arithmetic_intensity);
  profiler_stats vort_stats = finish_profiler(vort_profiler);
  printf("- Vort Calculation: %4.2f Flops/Cycle, %10ld cycles in %d runs. "
         "Arithmetic intensity: %4.2f\n",
         vort_stats.performance, vort_stats.cycles, vort_stats.runs, vort_stats.arithmetic_intensity);
#endif

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

#ifdef DEBUG
  char timestamp_filename[100];
  sprintf(timestamp_filename, "%s/timestamp_%d_%d_%d.txt", folder_name.c_str(), Nx, Ny, Nt);
  FILE *timestamp_file = fopen(timestamp_filename, "w");
  fclose(timestamp_file);
  make_latest_output(folder_name);
#endif

  return 0;
}

int main(int argc, char const *argv[]) {
  unsigned long long start_cycle, end_cycle;
  time_t start_sec, end_sec;
  asm volatile("RDTSC" : "=A"(start_cycle));
  time(&start_sec);
  run();
  asm volatile("RDTSC" : "=A"(end_cycle));
  time(&end_sec);
  printf("Cycles taken: %llu (%ld seconds)\n", end_cycle - start_cycle, end_sec - start_sec);

  return 0;
}
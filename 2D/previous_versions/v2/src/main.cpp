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
#include <stdexcept>
#include <string>
#include <vector>

#include "include/npy.hpp"
#include "include/profile.h"
#include "include/utils.h"

#ifdef BENCHMARK
#include <papi.h>
#endif

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

int Nx = 400;           // resolution in x
int Ny = 100;           // resolution in y
int Nt = 5000;          // number of timesteps
#define rho0 0.01       // reciprocal average density
#define tau -1.66666667 // reciprocal collision timescale (1/0.6)

#ifdef BENCHMARK
int papi_event_set = PAPI_NULL;

long long papi_values[2];

long long papi_drift_values[2] = {0, 0};
long long papi_rho_values[2] = {0, 0};
long long papi_feq_values[2] = {0, 0};
long long papi_f_values[2] = {0, 0};
long long papi_vort_values[2] = {0, 0};
#endif

// #define Nt 30  // number of timesteps

// Lattice speeds / weights
#define NL 9
const double idx[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
const double cxs[9] = {0.0, 0.0, 1.0, 1.0, 1.0, 0.0, -1.0, -1.0, -1.0};
const double cys[9] = {0.0, 1.0, 1.0, 0.0, -1.0, -1.0, -1.0, 0.0, 1.0};
const double weights[9] = {4.0 / 9, 1.0 / 9, 1.0 / 36, 1.0 / 9, 1.0 / 36, 1.0 / 9, 1.0 / 36, 1.0 / 9, 1.0 / 36}; // sums to 1

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

#ifdef BENCHMARK
  if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT) {
    fprintf(stderr, "PAPI library init error!\n");
    exit(1);
  }

  if (PAPI_create_eventset(&papi_event_set) != PAPI_OK) {
    fprintf(stderr, "PAPI create event set error!\n");
    exit(1);
  }

  // Register memory event
  if (PAPI_add_named_event(papi_event_set, "ANY_DATA_CACHE_FILLS_FROM_SYSTEM:MEM_IO_LCL") != PAPI_OK) {
    fprintf(stderr, "Couldn't add memory event to PAPI!\n");
    exit(1);
  }
  // Register flops event
  if (PAPI_add_event(papi_event_set, PAPI_FP_OPS) != PAPI_OK) {
    fprintf(stderr, "Couldn't add flop event to PAPI\n");
    exit(1);
  }
#endif

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

  double *Feq = (double *)malloc(Ny * Nx * NL * sizeof(double));
  double *F = (double *)malloc(Ny * Nx * NL * sizeof(double));
  double *vorticity = (double *)malloc(Ny * Nx * sizeof(double));
  double *rho = (double *)calloc(Ny * Nx, sizeof(double));
  double *cylinder = (double *)malloc(Ny * Nx * sizeof(double));
  double *ux = (double *)malloc(Ny * Nx * sizeof(double));
  double *uy = (double *)malloc(Ny * Nx * sizeof(double));
  double *temp = (double *)malloc(Ny * Nx * sizeof(double));
  int *x_coords = (int *)malloc(Ny * Nx * sizeof(int));
  int *y_coords = (int *)malloc(Ny * Nx * sizeof(int));

  debug_print("Initializing\n");

  srand(42); // some seed

  // Initialize F
  // flops = Ny*Nx*NL(1 div + 1 add + 1 add + 1 mult)
  for (int i = 0; i < Ny; i++) {
    for (int j = 0; j < Nx; j++) {
      for (int k = 0; k < NL; k++) {
        double rand_val = ((double)rand() / (RAND_MAX)) + 1;
        F[i * (Nx * NL) + j * NL + k] = 1 + 0.01 * rand_val;
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
      F[i * (Nx * NL) + j * NL + 1] += 2.0 * (1.0 + 0.2 * cos(2.0 * M_PI * (double)x_coords[i * Nx + j] / (double)Nx * 4.0));
    }
  }

  // flops = Ny*Nx*NL*(1 add)
  int res = 0;
  for (int i = 0; i < Ny; i++) {
    for (int j = 0; j < Nx; j++) {
      res = 0;
      for (int k = 0; k < NL; k++) {
        res += F[i * (Nx * NL) + j * NL + k];
      }
      rho[i * Nx + j] = res;
    }
  }

  // 	for i in idxs:		F[:,:,i] *= rho0 / rho
  // flops = Ny*Nx*NL*(2 adds)
  for (int j = 0; j < Ny; j++) {
    for (int k = 0; k < Nx; k++) {
      for (int i = 0; i < NL; i++) {
        F[j * (Nx * NL) + k * NL + i] *= rho0 * rho[j * Nx + k];
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

  double *bndryF = (double *)malloc(bndryF_size * NL * sizeof(double));

  // Currently assuming that every read/write is a miss
  profiler *rho_profiler = init_profiler(5 * Ny * Nx * NL + 3 * Ny * Nx, 8 * 5 * Ny * Nx * NL + 3 * Ny * Nx), *feq_profiler = init_profiler(20 * Nx * Ny * NL, 8 * 13 * Nx * Ny * NL), *f_profiler = init_profiler(3 * Nx * Ny * NL, 8 * 3 * Nx * Ny * NL),
           *vort_profiler = init_profiler(3 * Nx * Ny, 8 * 6 * Nx * Ny);

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
            temp[((k + shiftY + Ny) % Ny) * Nx + ((l + shiftX + Nx) % Nx)] = F[k * (Nx * NL) + l * NL + j];
          }
        }

        for (int k = 0; k < Ny; k++) {
          for (int l = 0; l < Nx; l++) {
            F[k * (Nx * NL) + l * NL + j] = temp[k * Nx + l];
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
              bndryF[index_bndryF * NL + l] = F[j * (Nx * NL) + k * NL + l];
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
#ifdef BENCHMARK
      if (PAPI_start(papi_event_set) != PAPI_OK) {
        fprintf(stderr, "PAPI start error!\n");
        exit(1);
      }
#endif
      start_run(rho_profiler);
      for (int j = 0; j < Ny; j++) {
        for (int k = 0; k < Nx; k++) {
          double res1 = 0;
          double res2 = 0;
          double res3 = 0;
          for (int l = 0; l < NL; l++) {
            res1 += F[j * (Nx * NL) + k * NL + l];
            res2 += F[j * (Nx * NL) + k * NL + l] * cxs[l];
            res3 += F[j * (Nx * NL) + k * NL + l] * cys[l];
          }
          double inv = 1 / res1;
          rho[j * Nx + k] = res1;
          ux[j * Nx + k] = res2 * inv;
          uy[j * Nx + k] = res3 * inv;
        }
      }
      end_run(rho_profiler);
#ifdef BENCHMARK
      if (PAPI_stop(papi_event_set, papi_values) != PAPI_OK) {
        fprintf(stderr, "PAPI stop error!\n");
        exit(1);
      }
      papi_rho_values[0] += 64 * papi_values[0];
      papi_rho_values[1] += papi_values[1];
#endif

      // set to zero
      memset(Feq, 0, Ny * Nx * NL * sizeof(double));

// flops = NL*Ny*Nx*(9 mults + 2 div  + 3 pows + 6 adds)
#ifdef BENCHMARK
      if (PAPI_start(papi_event_set) != PAPI_OK) {
        fprintf(stderr, "PAPI start error!\n");
        exit(1);
      }
#endif
      start_run(feq_profiler);
      for (int k = 0; k < NL; k++) {
        for (int j = 0; j < Ny; j++) {
          for (int l = 0; l < Nx; l++) {

            double rho_val = rho[j * Nx + l];
            double weight_val = weights[k];

            // 3*(cx*ux+cy*uy)
            double first = 3 * (cxs[k] * ux[j * Nx + l] + cys[k] * uy[j * Nx + l]);

            // 9*(cx*ux+cy*uy)**2/2
            double second = 9 * pow(cxs[k] * ux[j * Nx + l] + cys[k] * uy[j * Nx + l], 2) / 2;

            // 3*(ux**2+uy**2)/2
            double third = 3 * (pow(ux[j * Nx + l], 2) + pow(uy[j * Nx + l], 2)) / 2;

            Feq[j * (Nx * NL) + l * NL + k] = rho_val * weight_val * (1 + first + second - third);
          }
        }
      }
      end_run(feq_profiler);
#ifdef BENCHMARK
      if (PAPI_stop(papi_event_set, papi_values) != PAPI_OK) {
        fprintf(stderr, "PAPI stop error!\n");
        exit(1);
      }
      papi_feq_values[0] += 64 * papi_values[0];
      papi_feq_values[1] += papi_values[1];
#endif

// F += -(1.0/tau) * (F - Feq)
// flops = Ny*Nx*NL (2 add + 1 mult )
#ifdef BENCHMARK
      if (PAPI_start(papi_event_set) != PAPI_OK) {
        fprintf(stderr, "PAPI start error!\n");
        exit(1);
      }
#endif
      start_run(f_profiler);
      for (int j = 0; j < Ny; j++) {
        for (int k = 0; k < Nx; k++) {
          for (int l = 0; l < NL; l++) {
            F[j * (Nx * NL) + k * NL + l] += tau * (F[j * (Nx * NL) + k * NL + l] - Feq[j * (Nx * NL) + k * NL + l]);
          }
        }
      }
      end_run(f_profiler);
#ifdef BENCHMARK
      if (PAPI_stop(papi_event_set, papi_values) != PAPI_OK) {
        fprintf(stderr, "PAPI stop error!\n");
        exit(1);
      }
      papi_f_values[0] += 64 * papi_values[0];
      papi_f_values[1] += papi_values[1];
#endif

      // Apply boundary
      // F[cylinder,:] = bndryF
      // flops = 0
      int index_bndryF2 = 0;
      for (int j = 0; j < Ny; j++) {
        for (int k = 0; k < Nx; k++) {
          if (cylinder[j * Nx + k] == 1) {
            for (int l = 0; l < NL; l++) {
              F[j * (Nx * NL) + k * NL + l] = bndryF[index_bndryF2 * NL + l];
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
#ifdef BENCHMARK
      if (PAPI_start(papi_event_set) != PAPI_OK) {
        fprintf(stderr, "PAPI start error!\n");
        exit(1);
      }
#endif
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
#ifdef BENCHMARK
      if (PAPI_stop(papi_event_set, papi_values) != PAPI_OK) {
        fprintf(stderr, "PAPI stop error!\n");
        exit(1);
      }
      papi_vort_values[0] += 64 * papi_values[0];
      papi_vort_values[1] += papi_values[1];
#endif

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
  printf("- Rho  Calculation: %4.2f Flops/Cycle, %10llu cycles in %d runs. "
         "Arithmetic intensity: %4.2f\n",
         rho_stats.performance, rho_stats.cycles, rho_stats.runs, rho_stats.arithmetic_intensity);
  profiler_stats feq_stats = finish_profiler(feq_profiler);
  printf("- FEQ  Calculation: %4.2f Flops/Cycle, %10llu cycles in %d runs. "
         "Arithmetic intensity: %4.2f\n",
         feq_stats.performance, feq_stats.cycles, feq_stats.runs, feq_stats.arithmetic_intensity);
  profiler_stats f_stats = finish_profiler(f_profiler);
  printf("- F    Calculation: %4.2f Flops/Cycle, %10llu cycles in %d runs. "
         "Arithmetic intensity: %4.2f\n",
         f_stats.performance, f_stats.cycles, f_stats.runs, f_stats.arithmetic_intensity);
  profiler_stats vort_stats = finish_profiler(vort_profiler);
  printf("- Vort Calculation: %4.2f Flops/Cycle, %10llu cycles in %d runs. "
         "Arithmetic intensity: %4.2f\n",
         vort_stats.performance, vort_stats.cycles, vort_stats.runs, vort_stats.arithmetic_intensity);
#endif

#ifdef BENCHMARK
  printf("Drift Mem Transfer: %lld\n", papi_drift_values[0] / Nt);
  printf("Drift Floating point operations: %lld\n", papi_drift_values[1] / Nt);
  printf("Drift Arithmetic Intensity: %f\n", (double)papi_drift_values[1] / (double)papi_drift_values[0]);

  printf("Rho Mem Transfer: %lld\n", papi_rho_values[0] / Nt);
  printf("Rho Floating point operations: %lld\n", papi_rho_values[1] / Nt);
  printf("Rho Arithmetic Intensity: %f\n", (double)papi_rho_values[1] / (double)papi_rho_values[0]);

  printf("FEQ Mem Transfer: %lld\n", papi_feq_values[0] / Nt);
  printf("FEQ Floating point operations: %lld\n", papi_feq_values[1] / Nt);
  printf("FEQ Arithmetic Intensity: %f\n", (double)papi_feq_values[1] / (double)papi_feq_values[0]);

  printf("F Mem Transfer: %lld\n", papi_f_values[0] / Nt);
  printf("F Floating point operations: %lld\n", papi_f_values[1] / Nt);
  printf("F Arithmetic Intensity: %f\n", (double)papi_f_values[1] / (double)papi_f_values[0]);

  printf("Vort Mem Transfer: %lld\n", papi_vort_values[0] / Nt);
  printf("Vort Floating point operations: %lld\n", papi_vort_values[1] / Nt);
  printf("Vort Arithmetic Intensity: %f\n", (double)papi_vort_values[1] / (double)papi_vort_values[0]);

  // Cleanup PAPI
  PAPI_shutdown();
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
// free(BIG_CHUNGUS);
#ifdef DEBUG
  make_latest_output(folder_name);
#endif

  return 0;
}

int main(int argc, char const *argv[]) {
  if (argc == 4) {
    Nx = atoi(argv[1]);
    Ny = atoi(argv[2]);
    Nt = atoi(argv[3]);
  } else if (argc != 1) {
    printf("Usage: %s [Nx Ny Nt]\n", argv[0]);
    return 1;
  }

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

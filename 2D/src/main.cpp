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
#define tau_plus_1 (tau + 1)

// Lattice speeds / weights
#define NL 9

#define tau_4_9 (tau * 4.0 / 9)
#define tau_1_9 (tau * 1.0 / 9)
#define tau_1_36 (tau * 1.0 / 36)
void meshgrid(int *x_coords, int *y_coords) {
  for (int i = 0; i < Ny; ++i) {
    for (int j = 0; j < Nx; ++j) {
      x_coords[i * Nx + j] = j;
      y_coords[i * Nx + j] = i;
    }
  }
}

inline int scalar_index(int x, int l) { return (x * Nx) + l; }
inline int scalar_index(int y, int x, int l) { return y * Nx * NL + x * Nx + l; }

double *Feq;
double *F;
double *vorticity;
double *rho;
bool *cylinder;
int *collision_shape;
double *ux;
double *uy;
double *F_temp;
int *x_coords;
int *y_coords;
int bndryF_size = 0;
double *bndryF;
profiler *rho_profiler = init_profiler(5 * Ny * Nx * NL + 3 * Ny * Nx, 8 * 5 * Ny * Nx * NL + 3 * Ny * Nx);
profiler *feq_profiler = init_profiler(Ny * Nx * 67, 8 * 13 * Nx * Ny * NL);
profiler *f_profiler = init_profiler(2 * Nx * Ny * NL, 8 * 3 * Nx * Ny * NL);
profiler *vort_profiler = init_profiler(3 * Nx * Ny, 8 * 6 * Nx * Ny);
profiler *drift_profiler = init_profiler(0, Nx *Ny * 9 * 2 * 8);

string folder_name;

#ifdef BENCHMARK
int papi_event_set = PAPI_NULL;

long long papi_values[2];

long long papi_drift_values[2] = {0, 0};
long long papi_rho_values[2] = {0, 0};
long long papi_feq_values[2] = {0, 0};
long long papi_f_values[2] = {0, 0};
long long papi_vort_values[2] = {0, 0};
#endif

void initialize() {
  folder_name = make_output_folder(); // TODO delete empty folders

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
  Feq = (double *)aligned_alloc(32, Ny * Nx * NL * sizeof(double));
  F = (double *)aligned_alloc(32, Ny * Nx * NL * sizeof(double));
  vorticity = (double *)aligned_alloc(32, Ny * Nx * sizeof(double));
  rho = (double *)aligned_alloc(32, Ny * Nx * sizeof(double));
  cylinder = (bool *)aligned_alloc(32, Ny * Nx * sizeof(bool));
  ux = (double *)aligned_alloc(32, Ny * Nx * sizeof(double));
  uy = (double *)aligned_alloc(32, Ny * Nx * sizeof(double));
  F_temp = (double *)aligned_alloc(32, Ny * Nx * NL * sizeof(double));
  x_coords = (int *)aligned_alloc(32, Ny * Nx * sizeof(int));
  y_coords = (int *)aligned_alloc(32, Ny * Nx * sizeof(int));

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

  // TODO: cylinder is a constant. it should be an array of indeces what would
  // be 1, not an array of 0s and 1s
  for (int i = 0; i < Ny; i++) {
    for (int j = 0; j < Nx; j++) {
      cylinder[i * Nx + j] = (pow((double)x_coords[i * Nx + j] - (double)Nx / 4, 2) + pow((double)y_coords[i * Nx + j] - (double)Ny / 2, 2)) < pow(Ny / 4, 2);
    }
  }

  // loop trough cylinder and count the number of 1s
  // flops = 0 just reding the values
  for (int i = 0; i < Ny; i++) {
    for (int j = 0; j < Nx; j++) {
      if (cylinder[i * Nx + j] == 1) {
        bndryF_size++;
      }
    }
  }

  bndryF = (double *)aligned_alloc(32, bndryF_size * NL * sizeof(double));
}

void do_drift() {
  // # Drift
  // flops = 0
  memcpy(F_temp, F, Nx * Ny * NL * sizeof(double));

  for (int y = 0; y < Ny; y++) {
    for (int x = 0; x < Nx; x++) {
      F[scalar_index(y, 1, (x + 1) % Nx)] = F_temp[scalar_index(y, 1, x)];
      F[scalar_index((y + 1) % Ny, 2, (x + 1) % Nx)] = F_temp[scalar_index(y, 2, x)];
      F[scalar_index((y + 1) % Ny, 3, x)] = F_temp[scalar_index(y, 3, x)];
      F[scalar_index(((y + 1) % Ny), 4, ((x - 1 + Nx) % Nx))] = F_temp[scalar_index(y, 4, x)];
      F[scalar_index(((y + Ny) % Ny), 5, ((x - 1 + Nx) % Nx))] = F_temp[scalar_index(y, 5, x)];
      F[scalar_index(((y - 1 + Ny) % Ny), 6, ((x - 1 + Nx) % Nx))] = F_temp[scalar_index(y, 6, x)];
      F[scalar_index(((y - 1 + Ny) % Ny), 7, ((x + Nx) % Nx))] = F_temp[scalar_index(y, 7, x)];
      F[scalar_index(((y - 1 + Ny) % Ny), 8, ((x + 1) % Nx))] = F_temp[scalar_index(y, 8, x)];
    }
  }

  // bndryF = F[cylinder,:]
  // flops = 0
  int index_bndryF = 0;
  for (int y = 0; y < Ny; y++) {
    for (int x = 0; x < Nx; x++) {
      if (cylinder[scalar_index(y, x)] == 1) {
        for (int l = 0; l < NL; l++) {
          bndryF[index_bndryF * NL + l] = F[scalar_index(y, l, x)];
        }
        index_bndryF++;
      }
    }
  }

  // 0,1,2,3,4,5,6,7,8  INDEXES
  // reorder columns bndryF = bndryF[:,[0,5,6,7,8,1,2,3,4]]
  // flops = 0
  for (int j = 0; j < bndryF_size; j++) {
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
}

void do_rho() {
  // flops = Ny*Nx*21
  // bytes = Ny*Nx*8*2+Ny*Nx*8*9

  for (int j = 0; j < Ny; j++) {
    for (int k = 0; k < Nx; k++) {
      double res1 =
          F[scalar_index(j, 0, k)] + F[scalar_index(j, 1, k)] + F[scalar_index(j, 2, k)] + F[scalar_index(j, 3, k)] + F[scalar_index(j, 4, k)] + F[scalar_index(j, 5, k)] + F[scalar_index(j, 6, k)] + F[scalar_index(j, 7, k)] + F[scalar_index(j, 8, k)];
      double res2 = F[scalar_index(j, 2, k)] + F[scalar_index(j, 3, k)] + F[scalar_index(j, 4, k)] - F[scalar_index(j, 6, k)] - F[scalar_index(j, 7, k)] - F[scalar_index(j, 8, k)];
      double res3 = F[scalar_index(j, 1, k)] + F[scalar_index(j, 2, k)] - F[scalar_index(j, 4, k)] - F[scalar_index(j, 5, k)] - F[scalar_index(j, 6, k)] + F[scalar_index(j, 8, k)];

      double inv = 1 / res1;
      rho[scalar_index(j, k)] = res1;
      ux[scalar_index(j, k)] = res2 * inv;
      uy[scalar_index(j, k)] = res3 * inv;
    }
  }
}

void do_feq() {
  // flops = Ny*Nx*67
  // bytes =
  for (int y = 0; y < Ny; y++) {
    for (int x = 0; x < Nx; x++) {
      double third = 1 - 1.5 * (ux[y * Nx + x] * ux[y * Nx + x] + uy[y * Nx + x] * uy[y * Nx + x]);
      double weight_val;
      double curr = 0;
      double first = 0;
      double second = 0;
      Feq[scalar_index(y, 0, x)] = rho[scalar_index(y, x)] * tau_4_9 * third;

      Feq[scalar_index(y, 1, x)] = rho[scalar_index(y, x)] * tau_1_9 * (3 * uy[scalar_index(y, x)] + 4.5 * uy[scalar_index(y, x)] * uy[scalar_index(y, x)] + third);

      curr = ux[scalar_index(y, x)] + uy[scalar_index(y, x)];
      Feq[scalar_index(y, 2, x)] = rho[scalar_index(y, x)] * tau_1_36 * (3 * curr + 4.5 * curr * curr + third);

      Feq[scalar_index(y, 3, x)] = rho[scalar_index(y, x)] * tau_1_9 * (3 * ux[scalar_index(y, x)] + 4.5 * ux[scalar_index(y, x)] * ux[scalar_index(y, x)] + third);

      curr = ux[scalar_index(y, x)] - uy[scalar_index(y, x)];
      Feq[scalar_index(y, 4, x)] = rho[scalar_index(y, x)] * tau_1_36 * (3 * curr + 4.5 * curr * curr + third);

      Feq[scalar_index(y, 5, x)] = rho[scalar_index(y, x)] * tau_1_9 * (-3 * uy[scalar_index(y, x)] + 4.5 * uy[scalar_index(y, x)] * uy[scalar_index(y, x)] + third);

      curr = -ux[scalar_index(y, x)] - uy[scalar_index(y, x)];
      Feq[scalar_index(y, 6, x)] = rho[scalar_index(y, x)] * tau_1_36 * (3 * curr + 4.5 * curr * curr + third);

      Feq[scalar_index(y, 7, x)] = rho[scalar_index(y, x)] * tau_1_9 * (-3 * ux[scalar_index(y, x)] + 4.5 * ux[scalar_index(y, x)] * ux[scalar_index(y, x)] + third);

      curr = uy[scalar_index(y, x)] - ux[scalar_index(y, x)];
      Feq[scalar_index(y, 8, x)] = rho[scalar_index(y, x)] * tau_1_36 * (3 * curr + 4.5 * curr * curr + third);
    }
  }
}

void do_f() {
  for (int y = 0; y < Ny; y++) {
    for (int l = 0; l < NL; l++) {
      for (int x = 0; x < Nx; x++) {
        F[scalar_index(y, l, x)] = tau_plus_1 * F[scalar_index(y, l, x)] - Feq[scalar_index(y, l, x)];
      }
    }
  }
}

void do_vort() {

  int index_bndryF2 = 0;
  for (int y = 0; y < Ny; y++) {
    for (int x = 0; x < Nx; x++) {
      if (cylinder[scalar_index(y, x)] == 1) {
        ux[scalar_index(y, x)] = 0;
        uy[scalar_index(y, x)] = 0;
        for (int l = 0; l < NL; l++) {
          F[scalar_index(y, l, x)] = bndryF[index_bndryF2 * NL + l];
        }
        index_bndryF2++;
      }
    }
  }
  // flops = Ny*Nx*(3 adds)
  // Calculate j = 0 boundary
  // Calculate k = 0 boundary
  double ux_roll = ux[Nx - 1] - ux[1];
  double uy_roll = uy[(Ny - 1) * Nx] - uy[Nx];
  vorticity[0] = cylinder[0] == 1 ? 0 : ux_roll - uy_roll;

  // Calculate k = [1, Nx - 2]
  for (int x = 1; x < Nx - 1; x++) {
    if (cylinder[x] == 0) {
      double ux_roll = ux[x - 1] - ux[x + 1];
      double uy_roll = uy[scalar_index(Ny - 1, x)] - uy[scalar_index(1, x)];
      vorticity[x] = ux_roll - uy_roll;
    } else
      vorticity[x] = 0;
  }
  // Calculate k = Nx - 1 boundary
  ux_roll = ux[Nx - 2] - ux[0];
  uy_roll = uy[(Ny - 1) * Nx + Nx - 1] - uy[Nx + Nx - 1];
  vorticity[Nx - 1] = cylinder[Nx - 1] == 1 ? 0 : ux_roll - uy_roll;

  // Calculate j = [1, Ny - 2]
  for (int j = 1; j < Ny; j++) {
    // Calculate k = 0 boundary
    double ux_roll = ux[scalar_index(j, Nx - 1)] - ux[scalar_index(j, 1)];
    double uy_roll = uy[scalar_index(j - 1, 0)] - uy[scalar_index(j + 1, 0)];
    vorticity[scalar_index(j, 0)] = cylinder[scalar_index(j, 0)] == 1 ? 0 : ux_roll - uy_roll;

    // Calculate k = [1, Nx - 2]
    for (int k = 1; k < Nx - 1; k++) {
      if (cylinder[j * Nx + k] == 0) {
        double ux_roll = ux[scalar_index(j, k - 1)] - ux[j * Nx + k + 1];
        double uy_roll = uy[scalar_index(j - 1, k)] - uy[scalar_index(j + 1, k)];
        vorticity[scalar_index(j, k)] = ux_roll - uy_roll;
      } else
        vorticity[scalar_index(j, k)] = 0;
    }

    // Calculate k = Nx - 1 boundary
    if (cylinder[j * Nx + Nx - 1] == 0) {
      ux_roll = ux[scalar_index(j, Nx - 2)] - ux[scalar_index(j, 0)];
      uy_roll = uy[scalar_index(j - 1, Nx - 1)] - uy[scalar_index(j + 1, Nx - 1)];
      vorticity[scalar_index(j, Nx - 1)] = ux_roll - uy_roll;
    } else
      vorticity[scalar_index(j, Nx - 1)] = 0;
  }

  // Calculate j = Ny - 1 boundary
  // Calculate k = 0 boundary
  ux_roll = ux[scalar_index(Ny - 1, Nx - 1)] - ux[scalar_index(Ny - 1, 1)];
  uy_roll = uy[scalar_index(Ny - 2, 0)] - uy[0];
  vorticity[scalar_index(Ny - 1, 0)] = cylinder[scalar_index(Ny - 1, 0)] == 1 ? 0 : ux_roll - uy_roll;

  // Calculate k = [1, Nx - 2]
  for (int k = 1; k < Nx - 1; k++) {

    double ux_roll = ux[scalar_index(Ny - 1, k - 1)] - ux[scalar_index(Ny - 1, k + 1)];
    double uy_roll = uy[scalar_index(Ny - 2, k)] - uy[k];
    vorticity[scalar_index(Ny - 1, k)] = cylinder[scalar_index(Ny - 1, k)] == 1 ? 0 : ux_roll - uy_roll;
  }

  // Calculate k = Nx - 1 boundary
  ux_roll = ux[scalar_index(Ny - 1, Nx - 2)] - ux[scalar_index(Ny - 1, 0)];
  uy_roll = uy[scalar_index(Ny - 2, Nx - 1)] - uy[Nx - 1];
  vorticity[scalar_index(Ny - 1, Nx - 1)] = cylinder[scalar_index(Ny - 1, Nx - 1)] == 1 ? 0 : ux_roll - uy_roll;
}
void do_timestep() {
// ----------------- DRIFT -----------------
#ifdef BENCHMARK
  if (PAPI_start(papi_event_set) != PAPI_OK) {
    fprintf(stderr, "PAPI start error!\n");
    exit(1);
  }
#endif
  start_run(drift_profiler);
  do_drift();
  end_run(drift_profiler);
#ifdef BENCHMARK
  if (PAPI_stop(papi_event_set, papi_values) != PAPI_OK) {
    fprintf(stderr, "PAPI stop error!\n");
    exit(1);
  }
  papi_drift_values[0] += 64 * papi_values[0];
  papi_drift_values[1] += papi_values[1];
#endif

// ----------------- RHO -----------------
#ifdef BENCHMARK
  if (PAPI_start(papi_event_set) != PAPI_OK) {
    fprintf(stderr, "PAPI start error!\n");
    exit(1);
  }
#endif
  start_run(rho_profiler);
  do_rho();
  end_run(rho_profiler);
#ifdef BENCHMARK
  if (PAPI_stop(papi_event_set, papi_values) != PAPI_OK) {
    fprintf(stderr, "PAPI stop error!\n");
    exit(1);
  }
  papi_rho_values[0] += 64 * papi_values[0];
  papi_rho_values[1] += papi_values[1];
#endif

// ----------------- FEQ -----------------
#ifdef BENCHMARK
  if (PAPI_start(papi_event_set) != PAPI_OK) {
    fprintf(stderr, "PAPI start error!\n");
    exit(1);
  }
#endif
  start_run(feq_profiler);
  do_feq();
  end_run(feq_profiler);
#ifdef BENCHMARK
  if (PAPI_stop(papi_event_set, papi_values) != PAPI_OK) {
    fprintf(stderr, "PAPI stop error!\n");
    exit(1);
  }
  papi_feq_values[0] += 64 * papi_values[0];
  papi_feq_values[1] += papi_values[1];
#endif

// ----------------- F -----------------
#ifdef BENCHMARK
  if (PAPI_start(papi_event_set) != PAPI_OK) {
    fprintf(stderr, "PAPI start error!\n");
    exit(1);
  }
#endif
  start_run(f_profiler);
  do_f();
  end_run(f_profiler);
#ifdef BENCHMARK
  if (PAPI_stop(papi_event_set, papi_values) != PAPI_OK) {
    fprintf(stderr, "PAPI stop error!\n");
    exit(1);
  }
  papi_f_values[0] += 64 * papi_values[0];
  papi_f_values[1] += papi_values[1];
#endif

// ----------------- VORT -----------------
#ifdef BENCHMARK
  if (PAPI_start(papi_event_set) != PAPI_OK) {
    fprintf(stderr, "PAPI start error!\n");
    exit(1);
  }
#endif
  start_run(vort_profiler);
  do_vort();
  end_run(vort_profiler);
#ifdef BENCHMARK
  if (PAPI_stop(papi_event_set, papi_values) != PAPI_OK) {
    fprintf(stderr, "PAPI stop error!\n");
    exit(1);
  }
  papi_vort_values[0] += 64 * papi_values[0];
  papi_vort_values[1] += papi_values[1];
#endif
}

inline int run() {
  initialize();

#ifdef PROFILE
  const int PROFILE_RUNS = 5;
  const int PROFILE_DIGITS = floor(log10(PROFILE_RUNS)) + 1;
  printf("\rRun %-*d/%d done", PROFILE_DIGITS, 0, PROFILE_RUNS);
  fflush(stdout);
  for (int profile_counter = 0; profile_counter < PROFILE_RUNS; profile_counter++) {
#endif
    // Simulation loop
    for (int i = 0; i < Nt; i++) {
      do_timestep();
      debug_printf("\r%d", i);

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
  printf("\nProfiling results:\n");
  profiler_stats drift_stats = finish_profiler(drift_profiler);
  printf("- Drift  Calculation: %4.2f Flops/Cycle, %10ld cycles in %d runs. "
         "Arithmetic intensity: %4.2f\n",
         drift_stats.performance, drift_stats.cycles, drift_stats.runs, drift_stats.arithmetic_intensity);
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
  free(F_temp);

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

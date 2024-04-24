#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>

using namespace std;

const int Nx = 400;    // resolution in x
const int Ny = 400;    // resolution in y
const int rho0 = 100;  // average density
const int tau = 0.6;   // collision timescale
const int Nt = 4000;   // number of timesteps

void meshgrid(int x_coords[][Ny], int y_coords[][Ny]) {
  // Fill coordinate matrices
  for (int i = 0; i < Ny; ++i) {
    for (int j = 0; j < Nx; ++j) {
      x_coords[i][j] = j;
      y_coords[i][j] = i;
    }
  }
}

void roll_array(int i, int cx, int cy, double F[][Ny][Nx]) {
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


// TODO CHECK TYPES IF ALL CORRECT

int main() {
  // Lattice Boltzmann Simulation in 2D

  // Lattice speeds / weights
  int NL = 9;
  int idx[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  int cxs[9] = {0, 0, 1, 1, 1, 0, -1, -1, -1};
  int cys[9] = {0, 1, 1, 0, -1, -1, -1, 0, 1};
  double w[9] = {4 / 9, 1 / 9,  1 / 36, 1 / 9, 1 / 36,
                 1 / 9, 1 / 36, 1 / 9,  1 / 36};  // sums to 1

  // Initial conditions
  double F[Nx][Ny][NL];  // distribution function

  srand(42); // some seed

  // Initialize F
  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < NL; k++) {
        F[i][j][k] = 1 + 0.01 * (rand() % 100);  // TODO right?
      }
    }
  }

  int x_coords[Nx][Ny];
  int y_coords[Nx][Ny];

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
  int cylinderX[Nx][Ny];
  int cylinderY[Nx][Ny];

  meshgrid(cylinderX, cylinderY);

  // Simulation
  for (int i = 0; i < Nt; i++) {
    printf("Timestep %d\n", i);

    // Drift IDK about this one
    for (int j = 0; j < NL; j++) {
      //roll_array(j, cxs[j], cys[j], F);
    }
  }

  return 0;
}

// to run and compile run clear && gcc -o 2D/2D 2D/main.cpp -lm && 2D/2D
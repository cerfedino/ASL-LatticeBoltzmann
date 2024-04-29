#include <string>
#include <vector>
#include "./npy.hpp"
#include "./utils.h"

void save_npy_3d_double(double ***array, int x, int y, int z, std::string filename) {
  
  // convert array to vector  
  // this is needed because npy::write_npy expects a vector and cant deal with pointers to arrays somehow

  std::vector<double> vec;

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

void save_npy_2d_double(double **array, int x, int y, std::string filename) {

  // convert array to vector  
  // this is needed because npy::write_npy expects a vector and cant deal with pointers to arrays somehow

  std::vector<double> vec;

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

void save_npy_2d_int(int **array, int x, int y, std::string filename) {

  // convert array to vector
  // this is needed because npy::write_npy expects a vector and cant deal with pointers to arrays somehow
  // NOTE THIS FUNCTION SHOULD LIKELY NOT BE USED AS INT STORING SOMEHOW DOESNT WORK AND SAVES INCORRECTLY THE DATA

  std::vector<int> vec;
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
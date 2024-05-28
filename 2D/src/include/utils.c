#include "./npy.hpp"
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>


#ifdef MNx
#define Nx (MNx)
#endif
#ifdef MNy
#define Ny (MNy)
#endif
#ifdef MNt
#define Nt (MNt)
#endif

#ifdef DEBUG


#define debug_printf(fmt, ...) fprintf(stdout, fmt, __VA_ARGS__)
#define debug_print(fmt) fprintf(stdout, fmt)

void save_npy_3d_double(double *array, int x, int y, int z,
                        std::string filename) {
  // convert array to vector
  // this is needed because npy::write_npy expects a vector and cant deal with
  // pointers to arrays somehow

  std::vector<double> vec;

  for (int i = 0; i < x; i++) {
    for (int j = 0; j < y; j++) {
      for (int k = 0; k < z; k++) {
        vec.push_back(array[i * (y * z) + j * (z) + k]);
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

void save_npy_2d_double(double *array, int x, int y, std::string filename) {
  // convert array to vector
  // this is needed because npy::write_npy expects a vector and cant deal with
  // pointers to arrays somehow

  std::vector<double> vec;

  for (int i = 0; i < x; i++) {
    for (int j = 0; j < y; j++) {
      vec.push_back(array[i * y + j]);
    }
  }

  npy::npy_data<double> d;
  d.data = vec;
  d.shape = {(unsigned long)x, (unsigned long)y};
  d.fortran_order = false; // optional

  const std::string path{filename};
  npy::write_npy(path, d);
}

void save_npy_2d_int(int *array, int x, int y, std::string filename) {
  // convert array to vector
  // this is needed because npy::write_npy expects a vector and cant deal with
  // pointers to arrays somehow NOTE THIS FUNCTION SHOULD LIKELY NOT BE USED AS
  // INT STORING SOMEHOW DOESNT WORK AND SAVES INCORRECTLY THE DATA

  std::vector<int> vec;
  for (int i = 0; i < x; i++) {
    for (int j = 0; j < y; j++) {
      vec.push_back(array[i * y + j]);
    }
  }

  npy::npy_data<int> d;
  d.data = vec;
  d.shape = {(unsigned long)x, (unsigned long)y};
  d.fortran_order = false; // optional

  const std::string path{filename};
  npy::write_npy(path, d);
}

/// @brief This function prepares the output main output folder and creates a
/// new folder with the current date, time and parameters
/// @return
std::string make_output_folder() {
  // create new folder of format YYYY_MM_DD_HH_MM_SS
  time_t now = time(0);
  tm *ltm = localtime(&now);

  // its this ugly because i cant be bothered to do const char shenanigans
  char folder_name[100];
  sprintf(folder_name, "output/%02d_%02d_%02d_%02d_%02d_%02d_%d_%d_%d",
          1900 + ltm->tm_year, 1 + ltm->tm_mon, ltm->tm_mday, ltm->tm_hour,
          ltm->tm_min, ltm->tm_sec, Nx,Ny,Nt);

  if (mkdir(folder_name, 0755) == -1) {
    debug_printf("Could not create output folder with name %s\n", folder_name);
    /*throw runtime_error("Could not create output folder with name " +
                        string(folder_name));*/
  }

  return std::string(folder_name);
}

void make_latest_output(std::string folder) {
  char folder_name[220];
  sprintf(folder_name, "output/00_latest_%d_%d_%d", Nx,Ny,Nt);
  // check if output/latest folder exists in the current folder
  // if it does exist delete it
  if (access(folder_name, F_OK) == 0) {
    char delete_command[227];
    sprintf(delete_command, "rm -rf %s", folder_name);
    system((const char *)delete_command);
  }

  // copy over the latest output folder to output/latest
  char command[230];
  sprintf(command, "cp -r %s %s", folder.c_str(), folder_name);
  system((const char *)command);
}


#else
#endif
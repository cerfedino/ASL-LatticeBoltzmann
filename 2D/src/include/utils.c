#include "./npy.hpp"
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

#ifdef BENCHMARK
#include <papi.h>
#endif

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


#ifdef BENCHMARK

std::string papi_error_code_to_string(int code) {
  switch (code)
  {
    case 0:
      return "PAPI_OK (No error)";
    case -1:
      return "PAPI_EINVAL (Invalid argument)";
    case -2:
      return "PAPI_ENOMEM (Insufficient memory)";
    case -3:
      return "PAPI_ESYS (A System/C library call failed)";
    case -4:
      return "PAPI_ECMP (Not supported by component)";
    case -5:
      return "PAPI_ESBSTR (Backwards compatibility)";
    case -6:
      return "PAPI_EBUG (Internal error, please send mail to the developers)";
    case -7:
      return "PAPI_ENOEVNT (Event does not exist)";
    case -8:
      return "PAPI_ECNFLCT (Event exists, but cannot be counted due to counter resource limitations)";
    case -9:
      return "PAPI_ENOTRUN (EventSet is currently not running)";
    case -10:
      return "PAPI_EISRUN (EventSet is currently counting)";
    case -11:
      return "PAPI_ENOEVST (No such EventSet Available)";
    case -12:
      return "PAPI_ENOTPRESET (Event in argument is not a valid preset)";
    case -13:
      return "PAPI_ENOCNTR (Hardware does not support performance counters)";
    case -14:
      return "PAPI_EMISC (Unknown error code)";
    case -15:
      return "PAPI_EPERM (Permission level does not permit operation)";
    case -16:
      return "PAPI_ENOINIT (PAPI hasn't been initialized yet)";
    case -17:
      return "PAPI_ENOCMP (Component Index isn't set)";
    case -18:
      return "PAPI_ENOSUPP (Not supported)";
    case -19:
      return "PAPI_ENOIMPL (Not implemented)";
    case -20:
      return "PAPI_EBUF (Buffer size exceeded)";
    case -21:
      return "PAPI_EINVAL_DOM (EventSet domain is not supported for the operation)";
    case -22:
      return "PAPI_EATTR (Invalid or missing event attributes)";
    case -23:
      return "PAPI_ECOUNT (Too many events or attributes)";
    case -24:
      return "PAPI_ECOMBO (Bad combination of features)";
    case -25:
      return "PAPI_ECMP_DISABLED (Component containing event is disabled)";
    case -26:
      return "PAPI_EDELAY_INIT (Delayed initialization component)";
    case -27:
      return "PAPI_EMULPASS (Event exists, but cannot be counted due to multiple passes required by hardware)";
    case 28:
      return "PAPI_NUM_ERRORS (Number of error messages specified in this API)";

    default:
      return "Unknown error code: " + std::to_string(code);
  }
}

void papi_init(int *papi_event_set) {
 int papi_returned_version = PAPI_library_init(PAPI_VER_CURRENT);
  if (papi_returned_version != PAPI_VER_CURRENT) {
    fprintf(stderr, "PAPI library init error, we got version %d but we expected %d\n", papi_returned_version, PAPI_VER_CURRENT);
    exit(1);
  } else {
    fprintf(stdout, "PAPI library initialized successfully.\n");
  }

  int create_event_code = PAPI_create_eventset(papi_event_set);
  if (create_event_code != PAPI_OK) {
    fprintf(stderr, "PAPI create event set error with code %d %s\n", create_event_code, papi_error_code_to_string(create_event_code).c_str());
  } else {
    fprintf(stdout, "PAPI event set created successfully.\n");
  }

  // Register memory event
  int add_event_mem_io_lcl = PAPI_add_named_event(*papi_event_set, "ANY_DATA_CACHE_FILLS_FROM_SYSTEM:MEM_IO_LCL");
  if (add_event_mem_io_lcl != PAPI_OK) {
    fprintf(stderr, "PAPI add event MEM_IO_LCL error! Got code: %d %s\n", add_event_mem_io_lcl, papi_error_code_to_string(add_event_mem_io_lcl).c_str());

    // Register memory event (alternative if first fails)
    int add_event_mem_LCL_dram = PAPI_add_named_event(*papi_event_set, "DATA_CACHE_REFILLS_FROM_SYSTEM:LS_MABRESP_LCL_DRAM");
    if (add_event_mem_LCL_dram != PAPI_OK) {
      fprintf(stderr, "PAPI add event LS_MABRESP_LCL_DRAM error! Got code: %d %s\n", add_event_mem_LCL_dram, papi_error_code_to_string(add_event_mem_LCL_dram).c_str());
    } else {
      fprintf(stdout, "PAPI (BACKUP) event LS_MABRESP_LCL_DRAM added successfully.\n");
    }
  } else {
    fprintf(stdout, "PAPI event MEM_IO_LCL added successfully.\n");
  }

  // Register flops event
  int add_event_papi_fp_ops = PAPI_add_event(*papi_event_set, PAPI_FP_OPS);
  if (add_event_papi_fp_ops != PAPI_OK) {
    fprintf(stderr, "PAPI add event PAPI_FP_OPS error! Got code: %d %s\n", add_event_papi_fp_ops, papi_error_code_to_string(add_event_papi_fp_ops).c_str());
  } else {
    fprintf(stdout, "PAPI event PAPI_FP_OPS added successfully.\n");
  }
}

#endif

#include "utils.h"
#include <fstream>
#include <iostream>
// Now Linux only.
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#ifdef BENCHMARK
#include <papi.h>
#endif

#ifdef MNx
#define NX (MNx)
#endif
#ifdef MNy
#define NY (MNy)
#endif
#ifdef MNz
#define NZ (MNz)
#endif
#ifdef MNt
#define NT (MNt)
#endif

int scalar_index(int x, int y, int z) { return (z * NX * NY) + (y * NX) + x; }
int scalar_index(int x, int y, int z, int w) { return (x + y * NX + z * NX * NY + w * NX * NY * NZ); }

void output_array(double *array) {
  std::cout << "x,y,z value" << std::endl;
  for (int x = 0; x < NX; x++) {
    for (int y = 0; y < NY; y++) {
      for (int z = 0; z < NZ; z++) {
        std::cout << x << "," << y << "," << z << ": " << array[scalar_index(x, y, z)] << std::endl;
      }
    }
  }
}

void output_density() { output_array(density_field); }

void output_velocity() {
  int z_index = 0;
  for (int x = 0; x < NX; x++) {
    for (int y = 0; y < NY; y++) {
      std::cout << velocity_field[scalar_index(x, y, z_index)*3] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void output_indices_file() {
  std::ofstream output("output/indices.csv");
  output << "x,y,z" << '\n';
  for (int i = 0; i < NX; i++) {
    for (int j = 0; j < NY; j++) {
      for (int k = 0; k < NZ; k++) {
        output << i << "," << j << "," << k << '\n';
      }
    }
  }
  std::cout << std::endl;
  output.close();
}

void output_f_array(double *f_array, int z_index) {
  for (int i = 0; i < direction_size; i++) {
    std::cout << "ans(:,:," << z_index << "," << (i + 1) << ")" << '\n';
    std::cout << '\n';
    for (int x = 0; x < NX; x++) {
      for (int y = 0; y < NY; y++) {
        std::cout << f_array[scalar_index(x, y, z_index, i)] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << '\n';
  }
  std::cout << std::endl;
}

void output_lbm_data(std::string filename, bool header) {
  std::ofstream output_stream;
  output_stream.open(filename, std::ofstream::out | std::ofstream::app);
  if (header) {
    output_stream << "p,u_x,u_y,u_z" << '\n';
  }
  for (int x = 0; x < NX; x++) {
    for (int y = 0; y < NY; y++) {
      for (int z = 0; z < NZ; z++) {
        output_stream << density_field[scalar_index(x, y, z)] << "," << velocity_field[scalar_index(x, y, z)*3] << "," << velocity_field[scalar_index(x, y, z)*3+1] << "," << velocity_field[scalar_index(x, y, z)*3+2] << '\n';
      }
    }
  }
  output_stream.close();
}



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
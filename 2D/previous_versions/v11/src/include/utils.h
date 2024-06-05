#ifndef UTILS_H
#define UTILS_H

#include <string>

void save_npy_3d_double(double *array, int x, int y, int z,
                        std::string filename);
void save_npy_2d_double(double *array, int x, int y, std::string filename);
void save_npy_2d_int(int *array, int x, int y, std::string filename);

std::string make_output_folder();
void make_latest_output(std::string folder);

void papi_init(int *papi_event_set);

#endif
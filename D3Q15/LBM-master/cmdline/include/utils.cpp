#include "utils.h"
#include <iostream>
#include <fstream>
//Now Linux only.
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

//RapidJSON files.
#include "document.h"
#include "writer.h"
#include "stringbuffer.h"
#include "csv.h"

int NX=0, NY=2, NZ=3;

double c_s;
double nu;
double tau;

double *density_field;
vector_3_double *velocity_field;
double *previous_particle_distributions;
double *particle_distributions;
int direction_size = 15;

inline bool file_exists (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}
inline int scalar_index(int x, int y, int z) {
  return (z * NX * NY) + (y * NX) + x;
}
inline int scalar_index(int x, int y, int z, int w) {
  return (x + y * NX + z * NX * NY + w * NX * NY * NZ);
}

int setup(int& save_every, double& c_s, double& tau, int& n_steps, double& gamma_dot){

    std::cout<<"woop"<<std::endl;

    if(!file_exists("options.json")) {
    std::cout << "Please ensure that options.json exists. If not, it can be obtained from root directory of GitHub repo." << '\n';
    return -1;
  	}
	
	//std::cout << "Do you want to clean the previous run? (1 - Yes, 0 - No): ";
	int choice = 1;
	if(choice == 1) {
		system("rm -rf output");
		system("mkdir output");
	}

	std::ifstream t("options.json");
	std::string str((std::istreambuf_iterator<char>(t)),
	                 std::istreambuf_iterator<char>());
	rapidjson::Document d;
	d.Parse(str.c_str());
	rapidjson::Value& save_every_value = d["save_every"];
	save_every = save_every_value.GetInt();
	std::cout << "Save every: " << save_every << '\n';
	auto grid_size = d["grid_size"].GetArray();
	NX = grid_size[0].GetInt();
	NY = grid_size[1].GetInt();
	NZ = grid_size[2].GetInt();
	std::cout << "Grid size: " << NX << "x" << NY << "x" << NZ << '\n';
	rapidjson::Value& m_c_s = d["c_s"];
	c_s = m_c_s.GetDouble();
	std::cout << "c_s (Speed of sound): " << c_s <<" "<< '\n';
	rapidjson::Value& m_tau = d["tau"];
	tau = m_tau.GetDouble();
	std::cout << "Tau: " << tau << '\n';
	if(tau < 0.5) {
	    //Section 3.5.5.1 with delta t=1.
	    std::cout << "Error: Tau must be greater than 0.5 for numerical stability." << '\n';
	    return -1;
	}
	
	    rapidjson::Value& m_gamma_dot = d["gamma_dot"];
	    gamma_dot = m_gamma_dot.GetDouble();
	std::cout << "Shear rate (gamma_dot): " << gamma_dot << '\n';
    rapidjson::Value& m_n_steps = d["n_steps"];
    n_steps = m_n_steps.GetInt();
	return 1;
}

void output_array(double *array) {
  std::cout << "x,y,z value" << std::endl;
  for (int x = 0; x < NX; x++) {
    for (int y = 0; y < NY; y++) {
      for (int z = 0; z < NZ; z++) {
        std::cout << x << "," << y << "," << z << ": "
                  << array[scalar_index(x, y, z)] << std::endl;
      }
    }
  }
}

void output_density() { output_array(density_field); }



void output_velocity() {
  int z_index = 0;
  for (int x = 0; x < NX; x++) {
    for (int y = 0; y < NY; y++) {
      std::cout << velocity_field[scalar_index(x, y, z_index)].x << " ";
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
        output_stream << density_field[scalar_index(x, y, z)] << ","
                      << velocity_field[scalar_index(x, y, z)].x << ","
                      << velocity_field[scalar_index(x, y, z)].y << ","
                      << velocity_field[scalar_index(x, y, z)].z << '\n';
      }
    }
  }
  output_stream.close();
}


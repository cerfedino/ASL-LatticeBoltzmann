#include <stdio.h>
#include <iostream>
#include <fstream>
#include <streambuf>
//Now Linux only.
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "LBM.hpp"
#include "vector3.hpp"
#include "include/utils.h"
//RapidJSON files.
#include "document.h"
#include "writer.h"
#include "stringbuffer.h"
#include "csv.h"

using namespace rapidjson;


inline bool file_exists (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

int main(int argc, char** argv) {
  	
	int save_every, NX, NY, NZ, n_steps;
	double gamma_dot, c_s,tau;
	std::string velocity_set, boundary_conditions;
	int success = setup(save_every, NX, NY, NZ, c_s, tau, n_steps, gamma_dot, velocity_set, boundary_conditions);
	std::cout<<save_every<<NX<<NY<<NZ<<std::endl;
	if(success != 1)
	{
		std::cout<<"Setup failed!"<<std::endl;
		return 0;
	}

	unsigned long long start_cycle, end_cycle;
  	time_t start_sec, end_sec;
  	asm volatile("RDTSC" : "=A"(start_cycle));
  	time(&start_sec);

	//setup LBM
	LBM *solver = new LBM(NX,NY,NZ, velocity_set, c_s, tau, boundary_conditions, gamma_dot);
	for(int i = 0; i < argc; i++) {
		if(std::string(argv[i]) == "generate_ic") {
			solver->output_lbm_data("ic.csv", true);
			std::cout << "Generated ic.csv" << '\n';
			return 0;
		}
	}

	if(file_exists("ic.csv")) {
		std::cout << "Loading initial conditions" << '\n';
		io::CSVReader<4> in("ic.csv");
		in.read_header(io::ignore_extra_column, "p","u_x","u_y","u_z");
		double density,u_x,u_y,u_z;
		for(int i = 0; i < NX; i++) {
			for(int j = 0; j < NY; j++) {
				for(int k = 0; k < NZ; k++) {
					in.read_row(density,u_x,u_y,u_z);
					solver->set_density(i,j,k,density);
					solver->set_velocity(i,j,k,u_x,u_y,u_z);
				}
			}
		}
		std::cout << "Loaded initial conditions" << '\n';
	} else {
		std::cout << "Using default of p=1 for all x,y,z and u(x,t=0)=0 for all x,y,z. (Steady state)" << '\n';
		std::cout << "If you wish to use your own initial conditions, please run the program but with command: generate_ic as a argument which will output ic.csv in format of p,u_x,u_y,u_z, assume indexes are incrementing i,j,k for i<NX,j<NY and k<NZ" << '\n';
	}

	//Equation 3.5 with delta t = 1, LBM Principles and Practice book.
	double viscosity = c_s * c_s * (tau - 0.5);
	std::cout << "Kinematic shear viscosity: " << viscosity << '\n';
    //Equation 4.4.49
    std::cout << "For velocity set D2Q9,D3Q15 and D3Q27, |u_max|<0.577\n";
	solver->output_lbm_data("output/0.csv");
	solver->output_indices_file();
	solver->output_velocity();
	int scale = 1;
	int runs = n_steps * scale * scale * scale;

	//start simulation
	for(int i = 0; i < runs; i = i + 1) {
		solver->perform_timestep();
		if((i+1) % save_every == 0) {
            double percentage = (double) (i + 1) / (double) (runs) * 100.0;
            std::cout << "Saving data - " << (i + 1) << "/" << runs << " (" << percentage << "%)" << '\n';
            solver->output_lbm_data("output/" + std::to_string(i + 1) + ".csv");
            solver->output_velocity(); 
        }
	}
	std::cout << std::endl;
	delete solver;
	asm volatile("RDTSC" : "=A"(end_cycle));
  	time(&end_sec);
	printf("Cycles taken: %llu (%ld seconds)\n", end_cycle - start_cycle, end_sec - start_sec);

	return 0;
}

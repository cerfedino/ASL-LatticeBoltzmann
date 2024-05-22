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
#include "csv.h"

inline bool file_exists (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}


#define save_every 20
#define c_s  0.577350269
#define velocity_set "D3Q15"
#define nu 0.16666666666
#define tau 1
#define gamma_dot 0.01
#define boundary_conditions "lees_edwards"

using namespace std;

int main(int argc, char const *argv[]) {
    if (argc != 5) {
        cout << "Usage: ./lbm <Nx> <Ny> <Nz> <Nt>" << '\n';
        return 1;
    }

    int Nx = std::stoi(argv[1]);
    int Ny = std::stoi(argv[2]);
    int Nz = std::stoi(argv[3]);
    int Nt = std::stoi(argv[4]);

    cout << "Nx: " << Nx << " Ny: " << Ny << " Nz: " << Nz << " Nt: " << Nt << endl;


    
    system("rm -rf output");
    system("mkdir output");


	std::cout << "Tau: " << tau << '\n';
	if(tau < 0.5) {
	    //Section 3.5.5.1 with delta t=1.
	    std::cout << "Error: Tau must be greater than 0.5 for numerical stability." << '\n';
	    return -1;
	}

    if(velocity_set != "D3Q15" && velocity_set != "D3Q27" && velocity_set != "D2Q9") {
		std::cout << "Error: Please specify a valid velocity set such as D3Q15,D3Q27 or D2Q9." << '\n';
		return -1;
	}

	if(boundary_conditions != "lees_edwards" && boundary_conditions != "periodic" && boundary_conditions != "couette") {
	    std::cout << "Errors: boundary_conditions in options.json can either be: periodic, Couette (D2Q9 only) or lees_edwards (Lees-Edwards Shear, Please see research paper by Alexander Wagner)";
	    return -1;
	}

	if(Nz != 1 && velocity_set == "D2Q9") {
	    std::cout << "Warning: Nz=1 for D2Q9.";
	    return -1;
	}

	LBM *solver = new LBM(Nx,Ny,Nz, velocity_set, c_s, tau, boundary_conditions, gamma_dot);
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
		for(int i = 0; i < Nx; i++) {
			for(int j = 0; j < Ny; j++) {
				for(int k = 0; k < Nz; k++) {
					in.read_row(density,u_x,u_y,u_z);
					solver->set_density(i,j,k,density);
					solver->set_velocity(i,j,k,u_x,u_y,u_z);
				}
			}
		}
		std::cout << "Loaded initial conditions" << '\n';
	} else {
		std::cout << "Using default of p=1 for all x,y,z and u(x,t=0)=0 for all x,y,z. (Steady state)" << '\n';
		std::cout << "If you wish to use your own initial conditions, please run the program but with command: generate_ic as a argument which will output ic.csv in format of p,u_x,u_y,u_z, assume indexes are incrementing i,j,k for i<Nx,j<Ny and k<Nz" << '\n';
	}

	//Equation 3.5 with delta t = 1, LBM Principles and Practice book.
	double viscosity = c_s * c_s * (tau - 0.5);
	std::cout << "Kinematic shear viscosity: " << viscosity << '\n';
    //Equation 4.4.49
    std::cout << "For velocity set D2Q9,D3Q15 and D3Q27, |u_max|<0.577\n";
	solver->output_lbm_data("output/0.csv");
	solver->output_indices_file();
	int scale = 1;
	int runs = Nt * scale * scale * scale;
	for(int i = 0; i < runs; i = i + 1) {
		solver->perform_timestep();
		if((i+1) % save_every == 0) {
            double percentage = (double) (i + 1) / (double) (runs) * 100.0;
            std::cout << "Saving data - " << (i + 1) << "/" << runs << " (" << percentage << "%)" << '\n';
            solver->output_lbm_data("output/" + std::to_string(i + 1) + ".csv");
            //solver->output_velocity();
        }
	}
    
    std::ofstream ofs("output/timestamp.txt");
    ofs.close();

	std::cout << std::endl;
	delete solver;
	return 0;
}

#ifndef __Structures_h__
#define __Structures_h__
#include "basicFunctions.h"


//Declaration of constants
const std::vector<double> EPSILON_NULL;
const long PERMISSION = 500000;
const int LANCZOS_SIZE = 5;
const uInt BAND_LANCZOS_MAX_ITERATIONS = 1500;
const uInt NUM_THREADS_USED = 0.9 * omp_get_max_threads();


struct hubbardParam{ //Hubbard parameters
	int n_sites;
	float u, mu;
	std::vector<double> t_matrix;

	//Used only for diagExact
	int N_e, S_z;

	void print() {
		cct("--- HUBBARD PARAMETERS ---",33);
		std::cout << "\nSites = " << n_sites << "\n";
		std::cout << "Hubbard block used : N" << N_e << "Sz" << S_z << "\n";
		std::cout << "U = " << u << "\tmu = " << mu << "\n";
		std::cout << "tMatrix = " << "\n";
		print_matrix(t_matrix.data(), n_sites, n_sites);
		std::cout << std::endl;
	}
};



struct Electrons{
	uShort up = 0, down = 0;
};

struct greenParam{ //Green parameters
	int g_added_site, g_added_spin;
	float g_eta;
	bool g_add_all, g_compute, g_average;

	void print(){
		cct("--- GREEN COMPUTE PARAMETERS ---",33);
		std::string truth;
		if (g_compute) truth = "TRUE";
		else truth = "FALSE";
		std::cout << "\nCompute? = " << truth <<std::endl;

		if (g_add_all) truth = "TRUE";
		else truth = "FALSE";
		std::cout << "Add all = " << truth << std::endl;
		if (!g_add_all) std::cout << "Added Site = " << g_added_site << "\n";
		
		std::cout << "Added Spin = " << g_added_spin << std::endl;
		std::cout << "Eta value = " << g_eta << std::endl;
		if (g_average) truth = "TRUE";
		else truth = "FALSE";
		std::cout << "Average = " << truth << "\n" << std::endl;
	}
};

struct samplingParam{ //Sampling parameters
	unsigned long reticle, sampling_size;
	
	//MH sampling
	float beta_MH;

	//H application for green
	uInt nHapply;
	float beta_Happly;

	//Fundamental truncation cutoff value
	float fund_tc = 1;
	
	std::vector<sType> init_state;

	void print() {
		cct("--- SAMPLING PARAMETERS ---",33);
		std::cout << "\nReticle = " << reticle << std::endl;
		std::cout << "Truncated Cutoff = " << fund_tc << std::endl;
		std::cout << "Sampling size: "<< sampling_size << std::endl;
		std::cout << "Beta_MH = " << beta_MH << std::endl;
		std::cout << "# H applied = " << nHapply << "\tBeta_Happly = " 
            << beta_Happly << std::endl;
		std::cout << "Initial states:" << std::endl;
		print_vector(init_state.data(), init_state.size());
		std::cout << std::endl;
	}
};

struct justManyVariables {
	struct hubbardParam hubP;
	struct samplingParam sP;
	struct greenParam gP;

};



#endif
 

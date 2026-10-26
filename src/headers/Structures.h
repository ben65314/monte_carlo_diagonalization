#ifndef __Structures_h__
#define __Structures_h__
#include "basicFunctions.h"


//Declaration of constants
const std::vector<double> EPSILON_NULL;
const long PERMISSION = 500;
const int LANCZOS_SIZE = 10;
const uInt BAND_LANCZOS_MAX_ITERATIONS = 1500;
//const uInt NUM_THREADS_USED = 1;
//
#ifdef _OPENMP
    const uInt NUM_THREADS_USED = 1 * omp_get_max_threads();
#else
    const uInt NUM_THREADS_USED = 1;
#endif // _OPENMP


struct hubbardParam{ //Hubbard parameters
	int n_sites;
	double u, mu;
	std::vector<double> t_matrix;

    std::vector<double> R;
    std::vector<double> K;
    std::vector<int> dimension = {0,0,0};
    sType DIM = 3;


	//Used only for diagExact
	int N_e, S_z;

	void print() {
		cct("--- HUBBARD PARAMETERS ---",33);
		std::cout << "\nSites = " << n_sites << "\n";
        std::cout << "Lattice ("<<dimension.at(0)<<"x"<<dimension.at(1)<<"x"<<dimension.at(2)<<")"<<std::endl;
        std::cout<<"R:"<<std::endl;
        print_matrix(R.data(),n_sites,(int)DIM,1,0);

		std::cout << "Hubbard block used : N" << N_e << "Sz" << S_z << "\n";
		std::cout << "U = " << u << "\tmu = " << mu << "\n";
		std::cout << "t-Matrix = " << "\n";
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
	uInt nHapply = 0;
	float beta_Happly = 0;

	//Fundamental truncation cutoff value
	float fund_tc = 1;

	std::vector<sType> init_state;

    //Quick filling modules
    /*
    tree_like_sampling will speed up the sampling process by never
    proposing a state that already was accepted. There is no garanty
    that the sample size will be met.
    */
    bool tree_like_sampling = false;
    /*
    pre_filling_layers will fill the bottom layers if they could be
    filled completely according to the given subspace sample size
    */
    bool pre_filling_layers = false;
    /*
    propose_filled_layer_state will allow to go back to any filled layer.
    If false, the sampling will be able to go back in the most recent filled
    layer but not the others
    */
    bool propose_filled_layer_state = true;


	void print() {
		cct("--- SAMPLING PARAMETERS ---",33);
		std::cout << "\nReticle = " << reticle << std::endl;
		std::cout << "Truncated Cutoff = " << fund_tc << std::endl;
		std::cout << "Sampling size: "<< sampling_size << std::endl;
		std::cout << "Beta_MH = " << beta_MH << std::endl;
		std::cout << "# H applied = " << nHapply << "\tBeta_Happly = "
            << beta_Happly << std::endl;
        std::cout << "Sampling add-ons"<<std::endl;
        std::cout << "\ttree_like_sampling ("<<tree_like_sampling<<")"<<std::endl;
        std::cout << "\tpre_filling_layers ("<<pre_filling_layers<<")"<<std::endl;
        std::cout << "\tpropose_filled_layer_state ("<<propose_filled_layer_state<<")"<<std::endl;
		std::cout << "Initial states:" << std::endl;
		print_vector(init_state.data(), init_state.size(),1,0);
		std::cout << std::endl;
	}
};

struct justManyVariables {
	struct hubbardParam hubP;
	struct samplingParam sP;
	struct greenParam gP;

};



#endif


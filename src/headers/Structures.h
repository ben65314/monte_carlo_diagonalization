#ifndef __Structures_h__
#define __Structures_h__
#include "basicFunctions.h"


//Declaration of constants
const std::vector<double> EPSILON_NULL;
const long PERMISSION = 500000;
const int LANCZOS_SIZE = 5;
const uInt BAND_LANCZOS_MAX_ITERATIONS = 1500;
const uInt NUM_THREADS_USED = 0.9 * omp_get_max_threads();
const double ALPHA_D = 1;
const double BETA_D = 0;

const std::complex<double> ALPHA = 1;
const std::complex<double> BETA = 0;

const int INCXY = 1;

struct hubbardParam{ //Hubbard parameters
	int n_sites;
	float u, mu;
	//struct tInfo tI;
	std::vector<double> tMatrix;

	//Used only for diagExact
	int nElec, totSpin;

	void print() {
		cct("--- HUBBARD PARAMETERS ---",33);
		std::cout<<"\nSites = "<<n_sites<<std::endl;
		std::cout<<"Hubbard block used : N"<<nElec<<"Sz"<<totSpin<<std::endl;
		std::cout<<"U = "<<u<<"\tmu = "<<mu<<std::endl;
		std::cout<<"tMatrix = "<<std::endl;
		print_matrix(tMatrix.data(), n_sites, n_sites);
		std::cout<<std::endl;
	}
};



struct Electrons{
	uShort up = 0, down = 0;
};

struct greenParam{ //Green parameters
	int g_AddedSite, g_AddedSpin;
	float g_EtaValue;
	bool g_AddAll, g_compute, g_average;

	void print(){
		cct("--- GREEN COMPUTE PARAMETERS ---",33);
		std::string truth;
		if (g_compute) truth = "TRUE";
		else truth = "FALSE";
		std::cout<<"\nCompute? = "<< truth <<std::endl;

		if (g_AddAll) truth = "TRUE";
		else truth = "FALSE";
		std::cout<<"Add all = "<< truth <<std::endl;
		if (!g_AddAll)std::cout<<"Added Site = "<<g_AddedSite<<std::endl;
		
		std::cout<<"Added Spin = "<<g_AddedSpin<<std::endl;
		std::cout<<"Eta value = "<<g_EtaValue<<std::endl;
		if (g_average) truth = "TRUE";
		else truth = "FALSE";
		std::cout<<"Average = "<< truth <<"\n"<<std::endl;
	}
};

struct samplingParam{ //Sampling parameters
	unsigned long reticle, samplingSize;
	
	//MH sampling
	float beta_MH;

	//H application for green
	uInt nHapply;
	float beta_Happly;

	//Fundamental truncation cutoff value
	float fund_tc = 1;
	
	std::vector<sType> initState;

	void print() {
		cct("--- SAMPLING PARAMETERS ---",33);
		std::cout<<"\nReticle = "<<reticle<<std::endl;
		std::cout<<"Truncated Cutoff = "<<fund_tc<<std::endl;
		std::cout<<"Sampling size: "<<samplingSize<<std::endl;
		std::cout<<"Beta_MH = "<<beta_MH<<std::endl;
		std::cout<<"# H applied = "<<nHapply<<"\tBeta_Happly = "<<beta_Happly<<std::endl;
		std::cout<<"Initial states:"<<std::endl;
		print_vector(initState.data(), initState.size());
		std::cout<<std::endl;
	}
};

struct justManyVariables {
	struct hubbardParam hubP;
	struct samplingParam sP;
	struct greenParam gP;

};



#endif
 

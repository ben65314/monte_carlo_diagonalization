#pragma once
#include "../src/headers/LanczosSolver.h"
#include "../src/headers/greenFunctions.h"

bool t_lanczosAlgorithm(bool details) {
	srand(0);
	bool status1, status2, allStatus;

	hubbardParam hubP;
	hubP.n_sites = 2;
	std::vector<double> t{0,-1,
						 -1,0};
	hubP.t_matrix = t;
	hubP.u=8;
	hubP.mu = 4;
	//Tests value
	std::vector<sType> sArr_vec({5,6,9,10});
	StatesR_T<sType,double> sArr(&sArr_vec, hubP);
	sArr.electrons = find_number_of_electron(sArr.get_at(0), hubP.n_sites);
    std::vector<double> init(4,0);
    int deg = 1;

    LanczosSolver<double, decltype(sArr)> LS;

	double energy = LS.lanczos_algorithm(&init, &sArr, &deg);

	//Goal values
	float energy_ref = -8.4721359;
    float fund_ref[4] = {0.16246, 0.688191, 0.688191, 0.16246};

	status1 =    (abs(abs(init[0]) - fund_ref[0]) < 10e-4) \
			  && (abs(abs(init[1]) - fund_ref[1]) < 10e-4)\
			  && (abs(abs(init[2]) - fund_ref[2]) < 10e-4)\
			  && (abs(abs(init[3]) - fund_ref[3]) < 10e-4);
	status2 = (abs(energy - energy_ref) < 10e-6);

	allStatus = status1 && status2;

	cct("lanczos_algorithm", 33); std::cout << " ------------------------------------";
	if (details) {

		std::cout << "\nTest energy 2 sites ----- "; is_success(status2);
		std::cout << "Test fund state 2 sites ----- "; is_success(status1);
		std::cout << "------------------------------------ ";
        cct("lanczos_algorithm", 33);
	}
	is_success(allStatus);
	return (allStatus);
}
bool t_bandLanczosAlgorithm(bool details){
	bool status1 = false, status2 = false, status3 = false, allStatus = false;
	//Preparation of hubbardParam
	hubbardParam hubP;
	hubP.n_sites = 4;
	std::vector<double> t{0,-1,0,-1,
						 -1,0,-1,0,
						  0,-1,0,-1,
						 -1,0,-1,0};
	hubP.t_matrix = t;
	hubP.u=8;
	hubP.mu = 4;

    LanczosSolver<double, StatesR_T<sType,double>> LS;

	//vector len 36
	std::vector<double> array_vectors(180,0);
	initial_vector((sType)36,array_vectors.data());
	initial_vector((sType)36,array_vectors.data()+36);
	initial_vector((sType)36,array_vectors.data()+36*2);
	initial_vector((sType)36,array_vectors.data()+36*3);
	initial_vector((sType)36,array_vectors.data()+36*4);

	std::vector<sType> sArr36{51,53,54,57,58,60,83,85,86,89,90,92,99,101,102,105,106,108,147,149,150,153,154,156,163,165,166,169,170,172,195,197,198,201,202,204};
	StatesR_T<sType,double> s36(&sArr36, hubP);

	std::vector<double> energies36;
	std::vector<double> sub_vectors36;
	std::vector<double> product36;

	unsigned int maxIter36 = 36;
	energies36 = LS.band_lanczos_algorithm(&array_vectors,5,36,&s36,&maxIter36,&sub_vectors36,&product36);

	//Vector len6
	std::vector<double> vector(18,0);

	initial_vector((sType)6,vector.data());
	initial_vector((sType)6,vector.data()+6);
	initial_vector((sType)6,vector.data()+12);

	std::vector<sType> sArr6{243,245,246,249,250,252};
	StatesR_T<sType,double> s6(&sArr6, hubP);

	std::vector<double> energies6;
	std::vector<double> sub_vectors6;
	std::vector<double> product6;
	unsigned int maxIter6 = 20;
	energies6 = LS.band_lanczos_algorithm(&vector,3,6,&s6,&maxIter6,&sub_vectors6, &product6);

	//vector len 24
	std::vector<sType> sArr24{55,59,61,62,87,91,93,94,103,107,109,110,151,155,157,158,167,171,173,174,199,203,205,206};
	StatesR_T<sType,double> s24(&sArr24, hubP);

	//s36.showAllStatesInorder();

    std::vector<double> fund(36,0);
    int deg = 1;
	LS.lanczos_algorithm(&fund, &s36, &deg);
	std::vector<double> array_vectors24(96,0);
	for (unsigned int i = 0; i < 4; i++){
		excited_vector_projection(true, i, 0, fund.data(), &s36, &s24, array_vectors24.data() + i * 24);
	}

	std::vector<double> energies24;
	std::vector<double> sub_vectors24;
	std::vector<double> product24;
	unsigned int maxIter24 = 24;
	energies24 = LS.band_lanczos_algorithm(&array_vectors24,4,24,&s24,&maxIter24,&sub_vectors24,&product24);

	double refE1 = -10, refE2 = -17.320234939844382, refE3 = -14.3245553203367582;

	status1 = (energies6.at(0)- refE1) < 10e-8;
	status2 = (energies36.at(0) - refE2) < 10e-8;
	status3 = (energies24.at(0) - refE3) < 10e-8;


	allStatus = status1 && status2 && status3;
	cct("band_lanczos_algorithm", 33); std::cout << " -------------------------------";
	if (details) {

		std::cout << "\nTest 6x6 energy ----------------------- "; is_success(status1);
		std::cout << "REF: "<<refE1<<"\tComputed value: "<<energies6.at(0)<<std::endl;
		std::cout << "Test 36x36 energy --------------------- "; is_success(status2);
		std::cout << "REF: "<<refE2<<"\tComputed value: "<<energies36.at(0)<<std::endl;
		std::cout << "Test 24x24 energy (e.excitation)------- "; is_success(status3);
		std::cout << "REF: "<<refE3<<"\tComputed value: "<<energies24.at(0)<<std::endl;
		std::cout << "------------------------------- "; cct("band_lanczos_algorithm", 33);
	}
	is_success(allStatus);

	return allStatus;
}
bool t_fundEnergy(bool details) {
	std::cout << std::fixed;
	std::cout << std::setprecision(16);
	bool status1, status2, allStatus;
	hubbardParam hubP;

	hubP.n_sites = 4;
	std::vector<double> t{0, -1, 0, -1, -1, 0, -1, 0, 0, -1, 0, -1, -1, 0, -1, 0};
	hubP.t_matrix = t;
	hubP.u=8;
	hubP.mu = 0;

	std::vector<sType> a{51,53,54,57,58,60,83,85,86,89,90,92,99,101,102,105,106,108,147,149,150,153,154,156,163,165,166,169,170,172,195,197,198,201,202,204};
	StatesR_T<sType,double> states(&a, hubP);

	std::vector<sType> b{5,6,9,10};
	hubbardParam hubP2;
	hubP2.n_sites = 2;
	std::vector<double> t2{0, -1, -1, 0};
	hubP2.t_matrix = t2;
	hubP2.u=8;
	hubP2.mu = 0;
	StatesR_T<sType,double> s2(&b, hubP2);

    LanczosSolver<double, decltype(s2)> LS;

	//Tests value
	double fE1, fE2;
    std::vector<double> fS1(4,0),fS2(36,0);
    int deg = 1;
	fE1 = LS.fund_energy(&fS1, &s2, &deg);

	fE2 = LS.fund_energy(&fS2, &states, &deg);
	//Goal values
	double fE1_ref = -0.472135954999, fE2_ref = -1.32023495827;

	status1 = (abs(fE1 - fE1_ref) < pow(10, -5));
	status2 = (abs(fE2 - fE2_ref) < pow(10, -5));

	allStatus = status1 && status2;
	cct("fund_energy", 33); std::cout << " ------------------------------------------";
	if (details) {
		std::cout << "\nTest 4X4 ------- "; is_success(status1);
		std::cout << "REF: "<<fE1_ref<<"\tComputed value: "<<fE1<<std::endl;
		std::cout << "Test 36X36 ----- "; is_success(status2);
		std::cout << "REF: "<<fE2_ref<<"\tComputed value: "<<fE2<<std::endl;
		std::cout << "------------------------------------------ ";
        cct("fund_energy", 33);
	}
	is_success(allStatus);
	return (allStatus);
}

bool t_computeH_phi_r(bool details) {
	double phi_n_full[4] = {1.08, 6.11, -0.48, 1.19};
	double phi_n_half[2] = {3.24, 5.12};

	double* h_phi_n_full = new double[4]();
	double* h_phi_n_half = new double[2]();

    hubbardParam hubP2;
    hubP2.n_sites = 2;
    hubP2.u = 8;
    hubP2.mu = 4;
    std::vector<double> mt2{0,-1,-1,0};
    hubP2.t_matrix = mt2;

	std::vector<sType> temp_full = {5,6,9,10};
	std::vector<sType> temp_half = {5,6};
	StatesR_T<sType,double> sArr_full(&temp_full, hubP2);
	StatesR_T<sType,double> sArr_half(&temp_half, hubP2);

    sArr_full.H(h_phi_n_full, phi_n_full);
    sArr_half.H(h_phi_n_half, phi_n_half);

	//Ref values
	double ref_full[4] = {-5.63, -51.15, 1.57, -5.63};
	double ref_half[2] = {-5.12, -44.38};

	bool status1 = false, status2 = false;
	for (uInt i = 0; i < 4; i++) {
		if(abs(ref_full[i]-h_phi_n_full[i]) > 10e-3) break;
		status1 = true;
	}

	for (uInt i = 0; i < 2; i++) {
		if(abs(ref_half[i]-h_phi_n_half[i]) > 10e-3) break;
		status2 = true;
	}
	bool allStatus = status1 && status2;
	cct("computeH_phi_r", 33); std::cout << " ---------------------------------------";
	if (details) {
	std::cout << "\nTest 2 sites 4X4 ------- "; is_success(status1);
	std::cout << "Test 2 sites 2X2 ----- "; is_success(status2);
	std::cout << "--------------------------------------- "; cct("computeH_phi_r", 33);
	}
	is_success(allStatus);
	delete[] h_phi_n_full;

	delete[] h_phi_n_half;
	return (allStatus);
}


bool allTests_lanczosAlgorithm(bool details) {
	bool result = true;
	cct("~~~~~~~~~~~~~~~~~~~~~~~~~ TEST LANCZOS ALGORITHM ~~~~~~~~~~~~~~~~~~~~~~~~~", 44);
	std::cout << std::endl;
	//result &= t_calculate_ab(details);
	result &= t_lanczosAlgorithm(details);
	result &= t_bandLanczosAlgorithm(details);
	result &= t_fundEnergy(details);
	result &= t_computeH_phi_r(details);
	//result *= t_fundEnergy2(true);
	cct("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~", 44);
	is_success(result);


	return 0;
}

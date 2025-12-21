#pragma once
#include "../src/headers/basicFunctions.h"
int ONE_BF = 1;


bool t_get_number_of_blocks(bool details) {
	bool status1, status2, status3, allStatus;
	int const S2 = 2, S3 = 3, S4 = 4;
	int const REF2 = 9, REF3 = 16, REF4 = 25;

	status1 = (get_number_of_blocks(S2) == REF2);
	status2 = (get_number_of_blocks(S3) == REF3);
	status3 = (get_number_of_blocks(S4) == REF4);

	allStatus = status1 && status2 && status3;

	cct("getNumberOfBlocks", 33); std::cout << " ------------------------------------";
	if (details) {

		std::cout << "\nTest 2 sites ----- "; is_success(status1);
		std::cout << "Test 3 sites ----- "; is_success(status2);
		std::cout << "Test 4 sites ----- "; is_success(status3);
		std::cout << "------------------------------------ "; cct("getNumberOfBlocks", 33);
	}

	is_success(allStatus);
	return (allStatus);
}
bool t_find_min_of_array(bool details) {
	bool status1, status2, status3, allStatus;


	std::vector<double> array1{2.4, 6.2, 2.34, 3.478, 3.4782};
	std::vector<double> array2{1673, 1256, 234, 257.2, 236.22, 2574, 674, 845, 345, 1043};
	std::vector<double> array3{5.000234, 5.000102, 5.000235, 5.000521};
	double minValue1, minValue2, minValue3;

	minValue1 = 2.34;
	minValue2 = 234;
	minValue3 = 5.000102;

	status1 = (minValue1 == find_min_of_array(array1));
	status2 = (minValue2 == find_min_of_array(array2));
	status3 = (minValue3 == find_min_of_array(array3));
	allStatus = status1 && status2 && status3;

	cct("findMinOfArray", 33); std::cout << " ---------------------------------------";
	if (details) {

		std::cout << "\nTest #1 ----- "; is_success(status1);
		std::cout << "Test #2 ----- "; is_success(status2);
		std::cout << "Test #3 ----- "; is_success(status3);
		std::cout << "--------------------------------------- "; cct("findMinOfArray", 33);
	}

	is_success(allStatus);
	return allStatus;
}
bool t_initial_vector(bool details) {
	bool status1, status2, allStatus;

	//Tests value
	sType size_of_vector = 3;
	double* test = new double[size_of_vector];
	initial_vector(size_of_vector, test);

	//Goal values
	sType ref_size_of_vector = 3;
	double* ref = new double[ref_size_of_vector];
	ref[0] = (float)rand() / RAND_MAX;
	ref[1] = (float)rand() / RAND_MAX;
	ref[2] = (float)rand() / RAND_MAX;

	status1 = !((abs(ref[0]-test[0])<10e-4)&&(abs(ref[1]-test[1])<10e-4)&&(abs(ref[2]-test[2])<10e-4));
	float norm = dnrm2_(&ref_size_of_vector, test, &ONE_BF);
	status2 = (1 == norm);
	allStatus = status1 && status2;

	cct("init_vector", 33); std::cout << " ------------------------------------------";
	if (details) {

		std::cout << "\nTest random vector ----- "; is_success(status1);
		std::cout << "Test normalisation ----- "; is_success(status2);
		std::cout << "------------------------------------------ "; cct("init_vector", 33);
	}
	is_success(allStatus);
    delete[] test;
    delete[] ref;
	return (allStatus);
}
bool t_normalize(bool details) {
	bool status1, status2, status3, allStatus;

    sType size2 = 2;
    sType size3 = 3;
    sType size6 = 6;

	//Tests value
	double test1[2] = {3, 4};
	double test2[3] = {5.2, 1.23, 7.32};
	double test3[6] = {0.82, 0.24, 0.26, 0.52, 0.98, 0.24};

	normalize(test1, size2);
	normalize(test2, size3);
	normalize(test3, size6);

	//Goal values
	double ref1[2] = {0.6, 0.8};
	double ref2[3] = {0.5737707698, 0.1357188552, 0.807692699};
	double ref3[6] = {0.5677494215, 0.1661705624, 0.1800181092, 0.3600362185, 0.6785297964, 0.1661705624};
	double minus = -1;

	daxpy_(&size2, &minus, ref1, &ONE_BF, test1, &ONE_BF);
	daxpy_(&size3, &minus, ref2, &ONE_BF, test2, &ONE_BF);
	daxpy_(&size6, &minus, ref3, &ONE_BF, test3, &ONE_BF);
	for (int i = 0; i < 2; i++){test1[i] = remove_zeros(test1[i]);}
	for (int i = 0; i < 3; i++){test2[i] = remove_zeros(test2[i]);}
	for (int i = 0; i < 6; i++){test3[i] = remove_zeros(test3[i]);}

	status1 = (dnrm2_(&size2, test1, &ONE_BF) < 10e-4);
	status2 = (dnrm2_(&size3, test2, &ONE_BF) < 10e-4);
	status3 = (dnrm2_(&size6, test3, &ONE_BF) < 10e-4);

	allStatus = status1 && status2 && status3;

	cct("normalize", 33); std::cout << " --------------------------------------------";
	if (details) {

		std::cout << "\nTest 1 ----- "; is_success(status1);
		std::cout << "Test 2 ----- "; is_success(status2);
		std::cout << "Test 3 ----- "; is_success(status3);
		std::cout << "-------------------------------------------- "; cct("normalize", 33);
	}
	is_success(allStatus);
	return (allStatus);
}
bool t_comb(bool details) {
	bool status1, status2, status3, allStatus;

	unsigned long comb1[] = { 7,3 };
	unsigned long comb2[] = { 24,14 };
	unsigned long comb3[] = { 42,12 };

	unsigned long comb1Result = 35;
	unsigned long comb2Result = 1961256;
	unsigned long comb3Result = 11058116888;

	unsigned long combValue1 = comb(comb1[0], comb1[1]);
	unsigned long combValue2 = comb(comb2[0], comb2[1]);
	unsigned long combValue3 = comb(comb3[0], comb3[1]);

	status1 = (comb1Result == combValue1);
	status2 = (comb2Result == combValue2);
	status3 = (comb3Result == combValue3);


	allStatus = status1 && status2 && status3;

	cct("comb", 33); std::cout << " -------------------------------------------------";
	if (details) {

		std::cout << "\nTest 7C3 ------- "; is_success(status1);
		std::cout<<"REF:"<<comb1Result<<"\tVALUE:"<<combValue1<<std::endl;
		std::cout << "Test 24C14 ----- "; is_success(status2);
		std::cout<<"REF:"<<comb2Result<<"\tVALUE:"<<combValue2<<std::endl;
		std::cout << "Test 42C12 ----- "; is_success(status3);
		std::cout<<"REF:"<<comb3Result<<"\tVALUE:"<<combValue3<<std::endl;
		std::cout << "------------------------------------------------- "; cct("comb", 33);
	}
	is_success(allStatus);


	return (allStatus);
}
bool t_calculate_sd(bool details) {
	bool status1, status2, status3, allStatus;

	std::vector<double> x{0.707, 0.3161, -0.916, -0.3675, -0.1235};

	std::vector<double> y{21.7434, 22.7597, 23.8877, 20.1325, 20.7718,
		21.0535, 23.7322, 21.7174, 20.0374, 20.8133};

	std::vector<double> z{1000.1851, 1001.5942, 1000.3816, 1000.5525, 1000.7508, 1001.2872, 1000.0466, 1000.9675, 1000.0743, 1000.2533, 1001.3206, 1001.8898, 1000.7542, 1000.447, 1000.8059};

	double xstd = 0.558508;
	double ystd = 1.315892;
	double zstd = 0.545276;


	status1 = (abs(xstd - calculate_sd(x)) < pow(10, -4));
	status2 = (abs(ystd - calculate_sd(y)) < pow(10, -4));
	status3 = (abs(zstd - calculate_sd(z)) < pow(10, -4));
	allStatus = status1 && status2 && status3;

	cct("calculateSD", 33); std::cout << " ------------------------------------------";
	if (details) {

		std::cout << "\nTest 1 ----- "; is_success(status1);
		std::cout << "Test 2 ----- "; is_success(status2);
		std::cout << "Test 3 ----- "; is_success(status3);
		std::cout << "------------------------------------------ "; cct("calculateSD", 33);
	}
	is_success(allStatus);
	return (allStatus);
}
bool t_combination_all(bool details) {
	bool status1, status2, status3, allStatus;

    // Test 4 sites
    int up = 2, down = 2, sites = 4;
    int nu_1 = 0, nu_2 = 1, nu_3 = 2;

    std::vector<sType> gen_states_1;
    std::vector<sType> gen_states_2;
    std::vector<sType> gen_states_3;

    combination_all(up, down, sites, nu_1, &gen_states_1);
    combination_all(up, down, sites, nu_2, &gen_states_2);
    combination_all(up, down, sites, nu_3, &gen_states_3);


    std::vector<sType> ref_states_1 = {105,60,90,165,195,150};
    std::vector<sType> ref_states_2 = {172,166,99,108,57,54,201,86,202,101,53,169,163,197,149,147,198,156,89,106,83,92,58,154};
    std::vector<sType> ref_states_3 = {85,102,51,153,170,204};

    status1 = (gen_states_1.size() == ref_states_1.size());
    for (unsigned int i = 0; i < ref_states_1.size(); i++){
        auto it = std::find(gen_states_1.begin(),gen_states_1.end(), ref_states_1.at(i));
        if (it == gen_states_1.end()) {
            status1 = false;
            break;
        }
    }
    status2 = (gen_states_2.size() == ref_states_2.size());
    for (unsigned int i = 0; i < ref_states_2.size(); i++){
        auto it = std::find(gen_states_2.begin(),gen_states_2.end(), ref_states_2.at(i));
        if (it == gen_states_2.end()) {
            status2 = false;
            break;
        }
    }
    status3 = (gen_states_3.size() == ref_states_3.size());
    for (unsigned int i = 0; i < ref_states_3.size(); i++){
        auto it = std::find(gen_states_3.begin(),gen_states_3.end(), ref_states_3.at(i));
        if (it == gen_states_3.end()) {
            status3 = false;
            break;
        }
    }

	allStatus = status1 && status2 && status3;

	cct("combination_all", 33); std::cout << " --------------------------------------";
	if (details) {

		std::cout << "\nTest 1 ----- "; is_success(status1);
		std::cout << "Test 2 ----- "; is_success(status2);
		std::cout << "Test 3 ----- "; is_success(status3);
		std::cout << "-------------------------------------- "; cct("combination_all", 33);
	}
	is_success(allStatus);
	return (allStatus);
}

bool allTests_basicFunctions(bool details) {
	bool result = 1;

	cct("~~~~~~~~~~~~~~~~~~~~~~~~~~ TEST BASIC FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~", 44);
	std::cout << std::endl;
	result &= t_get_number_of_blocks(details);
	result &= t_find_min_of_array(details);
	result &= t_initial_vector(details);
	result &= t_normalize(details);
	result &= t_comb(details);
	result &= t_calculate_sd(details);
	result &= t_combination_all(details);
	cct("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~", 44);
	is_success(result);

	return result;
}

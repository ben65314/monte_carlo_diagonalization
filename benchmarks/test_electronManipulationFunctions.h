#pragma once
#include "../src/headers/electronManipulationFunctions.h"

//States related basic functions
bool t_Hu(bool details) {
	bool status1, status2, status3, allStatus;

	long A = 13, B = 27, C = 158;
	int ref1 = 1, ref2 = 2, ref3 = 1;

	status1 = (Hu(A,2) == ref1);
	status2 = (Hu(B,3) == ref2);
	status3 = (Hu(C,4) == ref3);

	allStatus = status1 && status2 && status3;

	cct("Hu", 33); std::cout << " ---------------------------------------------------";
	if (details) {

		std::cout << "\nTest 1 ----- "; is_success(status1);
		std::cout << "Test 2 ----- "; is_success(status2);
		std::cout << "Test 3 ----- "; is_success(status3);
		std::cout << "--------------------------------------------------- ";
        cct("Hu", 33);
	}
	is_success(allStatus);
	return (allStatus);
}
bool t_c(bool details) {
	bool status1, status2, status3, status4, allStatus;

	//States A(13, 1, 2,'S');
	sType A = 5;
	sType B = 27, C = 158;
	//States ref1(15, 1, 2,'S');
	sType ref1 = 6;
	sType ref2 = 19, ref3 = 222, ref4 = 214;

	c_dag_operator(&A,0);
	c_operator(&A,1);
	c_dag_operator(&B,3);
	c_operator(&C,6);
	status3 = ((C == ref3)&&(C == ref3));
	c_dag_operator(&C,3);


	status1 = ((A == ref1)&&(A == ref1));
	status2 = (B == ref2);
	status4 = ((C == ref4)&&(C == ref4));

	allStatus = status1 && status2 && status3 && status4;

	cct("C", 33); std::cout << " ----------------------------------------------------";
	if (details) {

		std::cout << "\nAdd Down -------- "; is_success(status1);
		std::cout << "Remove Up ------- "; is_success(status2);
		std::cout << "Add Up ---------- "; is_success(status3);
		std::cout << "Remove Down ----- "; is_success(status4);
		std::cout << "---------------------------------------------------- ";
        cct("C", 33);
	}
	is_success(allStatus);
	return (allStatus);
}
bool t_Ht(bool details) {
	bool status1 = true, status2 = true, status3 = true, allStatus;

	//Tests value
	sType s1 = 10, s2 = 22, s3 = 187;

	//Parameters
	hubbardParam hubP2,hubP3,hubP4;
	hubP2.n_sites=2;
	hubP3.n_sites=3;
	hubP4.n_sites=4;
	std::vector<double> mt2{0,-1,-1,0};
	std::vector<double> mt3{0,-1,-1,-1,0,-1,-1,-1,0};
	std::vector<double> mt4{0,-1,0,-1,-1,0,-1,0,0,-1,0,-1,-1,0,-1,0};

	hubP2.t_matrix = mt2;
	hubP3.t_matrix = mt3;
	hubP4.t_matrix = mt4;
	//Goal values
	std::vector<sType> nextS1{6, 9};
	std::vector<sType> nextS2{14, 19, 21, 38};
	std::vector<sType> nextS3{123, 183, 189, 219};

	std::vector<sType> s1_;
	Ht(s1, &s1_, &hubP2);

	std::vector<sType> s2_;
	Ht(s2, &s2_, &hubP3);

	std::vector<sType> s3_;
	Ht(s3, &s3_, &hubP4);
	std::sort(s1_.data(), s1_.data() + s1_.size());
	std::sort(s2_.data(), s2_.data() + s2_.size());
	std::sort(s3_.data(), s3_.data() + s3_.size());
	if (s1_.size() == nextS1.size()) {
		for (unsigned long i = 0; i < s1_.size(); i++)
		{
			if (s1_.at(i) != nextS1.at(i)) {
				status1 = false;
			}
		}
	}
	else {
		status1 = false;
	}

	if (s2_.size() == nextS2.size()) {
		for (unsigned long i = 0; i < s2_.size(); i++)
		{
			if (s2_.at(i) != nextS2.at(i)) {
				status2 = false;
			}
		}
	}
	else {
		status2 = false;
	}

	if (s3_.size() == nextS3.size()) {
		for (unsigned long i = 0; i < s3_.size(); i++)
		{
			if (s3_.at(i) != nextS3.at(i)) {
				status3 = false;
			}
		}
	}
	else {
		status3 = false;
	}

	allStatus = status1 && status2 && status3;
	cct("Ht", 33); std::cout << " ---------------------------------------------------";
	if (details) {
		std::cout << "\nTest 2 sites ----- "; is_success(status1);
		print_vector(s1_.data(),s1_.size());
		std::cout<<" REF:";print_vector(nextS1.data(),nextS1.size());
        std::cout<<std::endl;
		std::cout << "Test 3 sites ----- "; is_success(status2);
		print_vector(s2_.data(),s2_.size());
		std::cout<<" REF:";print_vector(nextS2.data(),nextS2.size());
        std::cout<<std::endl;
		std::cout << "Test 4 sites ----- "; is_success(status3);
		print_vector(s3_.data(),s3_.size());
		std::cout<<" REF:";print_vector(nextS3.data(),nextS3.size());
        std::cout<<std::endl;
		std::cout << "--------------------------------------------------- ";
        cct("Ht", 33);
	}
	is_success(allStatus);
	return (allStatus);
}
bool t_JumpEnergy(bool details) {
	bool status1 = true, status2 = true, status3 = true, allStatus;

	//Tests value
	//Left states
	sType s1_l = 10, s2_l = 59, s3_l = 187;
	//Right states
	sType s1_r = 9, s2_r = 62, s3_r = 219;

	//Parameters
	hubbardParam hub2;
	hub2.n_sites = 2;
	hub2.u = 8;
	hub2.mu = 4;
	hub2.t_matrix = {0,-1,-1,0};

	hubbardParam hub3;
	hub3.n_sites = 3;
	hub3.u = 8;
	hub3.mu = 4;
	hub3.t_matrix = {0,-1,-1,-1,0,-1,-1,-1,0};

	hubbardParam hub4;
	hub4.n_sites = 4;
	hub4.u = 8;
	hub4.mu = 4;
	hub4.t_matrix = {0,-1,0,-1,-1,0,-1,0,0,-1,0,-1,-1,0,-1,0};

	//Ref values
	double ref1 = -1;
	double ref2 = 1;
	double ref3 = -1;

	std::vector<sType> proj;
	std::vector<double> energies;
	t_jump_energy(s1_r, &proj, &energies, &hub2);
	double tJumpElement = energies.at(std::find(proj.begin(), proj.end(), s1_l)
                                   - proj.begin());

	proj.clear(); energies.clear();
	t_jump_energy(s2_r, &proj, &energies, &hub3);
	double tJumpElement2 = energies.at(std::find(proj.begin(), proj.end(), s2_l)
                                    - proj.begin());

	proj.clear(); energies.clear();
	t_jump_energy(s3_r, &proj, &energies, &hub4);
	double tJumpElement3 = energies.at(std::find(proj.begin(), proj.end(), s3_l)
                                    - proj.begin());

	status1 = (ref1 == tJumpElement);
	status2 = (ref2 == tJumpElement2);
	status3 = (ref3 == tJumpElement3);

	allStatus = status1 && status2 && status3;
	cct("tJumpEnergy", 33); std::cout << " ------------------------------------------";
	if (details) {
		std::cout << "\nTest 2 sites ----- "; is_success(status1);
		std::cout<<" REF:"<<ref1<<"\tVALUE:"<<tJumpElement<<std::endl;
		std::cout << "Test 3 sites ----- "; is_success(status2);
		std::cout<<" REF:"<<ref2<<"\tVALUE:"<<tJumpElement2<<std::endl;
		std::cout << "Test 4 sites ----- "; is_success(status3);
		std::cout<<" REF:"<<ref3<<"\tVALUE:"<<tJumpElement3<<std::endl;
		std::cout << "------------------------------------------ ";
        cct("tJumpEnergy", 33);
	}
	is_success(allStatus);
	return (allStatus);
}
bool t_compute_mu(bool details) { //TO DO
	bool status1, status2, status3, allStatus;

	hubbardParam hub3;
	hub3.n_sites = 3;
	hub3.u = 8;
	hub3.mu = 3;
	hub3.t_matrix = {0,-1,-1,-1,0,-1,-1,-1,0};

	Electrons elec_value_1 = find_number_of_electron(10, hub3.n_sites);
	Electrons elec_value_2 = find_number_of_electron(11, hub3.n_sites);
	Electrons elec_value_3 = find_number_of_electron(25, hub3.n_sites);

	double refMu1 = -6, refMu2 = -9, refMu3 = -9;

	status1 = (refMu1 == compute_mu(hub3.mu, elec_value_1));
	status2 = (refMu2 == compute_mu(hub3.mu, elec_value_2));
	status3 = (refMu3 == compute_mu(hub3.mu, elec_value_3));

	allStatus = status1 && status2 && status3;

	cct("compute_mu", 33); std::cout << " -------------------------------------------";
	if (details) {
		std::cout << "\nTest 1 ----- "; is_success(status1);
		std::cout << "Test 2 ----- "; is_success(status2);
		std::cout << "Test 3 ----- "; is_success(status2);
		std::cout << "------------------------------------------- ";
        cct("compute_mu", 33);
	}
	is_success(allStatus);
	return (allStatus);
}
bool t_create_anti_ferro(bool details) {
	bool status1 = true, status2 = true, status3 = true, allStatus;

	//Tests value
    int n_up_1 = 2, n_up_2 = 1, n_up_3 = 3;
    int n_down_1 = 2, n_down_2 = 2, n_down_3 = 2;
    //Ref
	sType s1_ref = 165;
	sType s2_ref = 37;
	sType s3_ref = 181;


	//Parameters
	hubbardParam hub4;
	hub4.n_sites = 4;
	hub4.u = 8;
	hub4.mu = 4;
	hub4.t_matrix = {0,-1,0,-1,-1,0,-1,0,0,-1,0,-1,-1,0,-1,0};

    //Tests
    sType val_1 = create_anti_ferro(hub4.n_sites, n_up_1, n_down_1);
    sType val_2 = create_anti_ferro(hub4.n_sites, n_up_2, n_down_2);
    sType val_3 = create_anti_ferro(hub4.n_sites, n_up_3, n_down_3);

    status1 = (s1_ref == val_1);
    status2 = (s2_ref == val_2);
    status3 = (s3_ref == val_3);

	allStatus = status1 && status2 && status3;
	cct("create_anti_ferro", 33); std::cout << " ------------------------------------";
	if (details) {
		std::cout << "\nTest 2 sites ----- "; is_success(status1);
		std::cout<<" REF:"<<s1_ref<<"\tVALUE:"<<val_1<<std::endl;
		std::cout << "Test 3 sites ----- "; is_success(status2);
		std::cout<<" REF:"<<s2_ref<<"\tVALUE:"<<val_2<<std::endl;
		std::cout << "Test 4 sites ----- "; is_success(status3);
		std::cout<<" REF:"<<s3_ref<<"\tVALUE:"<<val_3<<std::endl;
		std::cout << "------------------------------------ ";
        cct("create_anti_ferro", 33);
	}
	is_success(allStatus);
	return (allStatus);
}
bool allTests_electronManipulationFunctions(bool details) {
	bool result = 1;

	cct("~~~~~~~~~~~~~~~~~~ TEST ELECTRON MANIPULATION FUNCTIONS ~~~~~~~~~~~~~~~~~~", 44);
	std::cout << std::endl;
	result &= t_Hu(details);
	result &= t_c(details);
	result &= t_Ht(details);
	result &= t_JumpEnergy(details);
	result &= t_compute_mu(details);
	result &= t_create_anti_ferro(details);
	cct("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~", 44);
	is_success(result);

	return result;
}

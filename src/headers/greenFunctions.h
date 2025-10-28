#include "LanczosSolver.h"

#ifndef __greenFunctions_h__
#define __greenFunctions_h__

template <class StatesArrType> void c_subspace(
        StatesArrType* const sub_space, int location, int spin, bool creation, 
        StatesArrType* c_sub_space) {
	/****************************************************************
	* Applies a c operator on an entire space of states
	* 
	* Parameters
	* ----------
	* sub_space	: (StatesArrType*) States who will be applied the c operator
	* location	: (int) site created or destroyed
	* spin		: (int) value of the spin created or destroyed (1=up, 0=down)
	* creation  : (bool) Is c a creation or annihilation operator
	* c_sub_space : (StatesArrType*) States resulting of the c appliacation
	*
	* Returns
	* -------
	* NONE
	*****************************************************************/
	unsigned int sites = sub_space->sys_hubP.n_sites;
	unsigned long len = sub_space->get_length();
	
	int index = sites - location - 1 + spin * sites;
    //Creates
	if (creation) {
		for (unsigned long i = 0; i < len; i++) {
			sType temp = sub_space->get_at(i);
			if (c_operator(&temp, index)) { 
				c_sub_space->add(temp);
			}
		}
	}
    //Annihilates 
	else {
		for (unsigned long i = 0; i < len; i++) {
			sType temp = sub_space->get_at(i);
			if (c_dag_operator(&temp, index)) { 
				c_sub_space->add(temp);
			}
		}
	}
}

template <class StatesArrType> void green_space_projection(
        StatesArrType* const origin_sub_space, int spin, bool creation, 
        StatesArrType* proj_sub_space) {
	/****************************************************************
	* Applies all the possible c_mu to create all the states of the projected space
	* 
	* Parameter
	* ----------
	* origin_sub_space : (StatesArrType) States who will be applied 
    *                                    the c operators
	* spin			   : (int) value of the spin created 
    *                          or destroyed (1=up, 0=down)
	* creation		   : (bool) Is c a creation or annihilation operator
	* proj_sub_space   : (StatesArrType) States resulting of the c application
    *
	* Templates:
	* ----------
    * StatesArrType    : StatesR_T, StatesR_H
	*
	* Returns
	* -------
	* NONE
	*****************************************************************/
	if (verbose < -9) std::cout << "green_space_projection(...)" << std::endl;

	unsigned int sites = origin_sub_space->sys_hubP.n_sites;
	proj_sub_space->remove_all();

	for (unsigned int i = 0; i < sites; i++) {
		c_subspace(origin_sub_space, i, spin, creation, proj_sub_space);
	}
}

std::string write_q_system(double fundE, double eta, Electrons elec,
                          hubbardParam* hubP_to_write, float perc_states) {
	/****************************************************************
	* Writes the system parameters for the qMatrices.txt file
	* 
	* Parameters
	* ----------
	* fundE			: (double) Fundamental Energy
	* eta			: (double) Width of the Lorentz
	* elec			: (Electrons) Electrons of the bloc system 
	* hubP_to_write : (hubbardParam*) System parameters to write
	* perc_states	: (float) {ercentage of states taken
	*
	* Returns
	* -------
	* text			: (std::string) String of parameters
	*****************************************************************/
	std::string text = "";
	text += to_string_p(hubP_to_write->n_sites, 0) + "\t";
	text += to_string_p(hubP_to_write->u, 0) + "\t";
	text += to_string_p(hubP_to_write->mu, 0) + "\t";
	text += to_string_p((elec.up+elec.down), 0) + "\t";
	text += to_string_p((elec.up-elec.down), 0) + "\t";
	text += to_string_p(fundE, 16) + "\t";
	text += to_string_p(eta, 5) + "\t";
	text += to_string_p(perc_states, 2) + "\t";
	return text;
}

template <class T> std::string write_q_matrix_string(
        std::vector<T>* QM, std::vector<double>* eigen_val, int ns) {
	/****************************************************************
	* Writes the Q-matrices of the system 
	* 
	* Parameter
	* ----------
	* QM		: (std::vector<T>*) Q-matrix
	* eigen_val	: (std::vector<double>*) eigen energies
	* ns		: (int) number of sites of the system
	*
	* Templates
	* ---------
	* T				: float, double, std::complex<double>
	*
	* Returns
	* -------
	* text			: (std::string) String of the Q-matrix
	*****************************************************************/
	std::string text = "";
	for (uInt j = 0; j < eigen_val->size(); j++) {
		text += to_string_pq(eigen_val->at(j),5,12) + "\t";
		for (int i = 0; i < ns; i++) {
			text += to_string_pq(QM->at(i * eigen_val->size() + j), 4, 14) 
                    + "\t";
		}
		if (j < eigen_val->size()-1) text+="\n";
	}

	return text;
}



template <class T, class StatesArrType> void excited_vector_projection(
        bool create, unsigned int site, int spin, const T* initial_vector, 
        const StatesArrType* initial_states, StatesArrType* projected_states, 
        T* projected_vector) {
	/****************************************************************
	* Transform a vector countained in a space to the according vector 
    * in the excited space. 
	* 
	* Parameters
	* ----------
	* create			: (bool) Are we creating or annihilating an electron 
    *                            to get to the excited space
	* site				: (int) which site is the electron or hole created
	* spin				: (int) value of the spin created 
    *                           or destroyed (1=up, 0=down)
	* initial_vector	: (const T*) the initial vector in 
    *                                the initial_states space
	* initial_states	: (const StatesArr*) the initial state space
	* projected_states	: (StatesArr*) the projected state space
	* projected_vector	: (T*) the projected vector in 
    *                               the projected_states space
	* 
	* Templates:
	* ----------
	* T					: double, std::complex<double>
    * StatesArrType     : StatesR_T, StatesR_H
	*
	* Returns
	* -------
	* NONE
	*****************************************************************/
	sType current_state;
	int ns = initial_states->sys_hubP.n_sites; //number of sites
	int index = ns - site - 1 + spin * ns;
	
	for(long unsigned i = 0; i < initial_states->get_length(); i++){
		//Exciting the state
		current_state = initial_states->get_at(i);
		bool can_projected = false;
		sType one = 1UL << (ns * 2 - 1);
		char phase = 1;
		for (int l = 2 * ns - 1; l > index; l--) {
			if ((one & current_state) != 0) phase *= -1;
			one >>= 1;
		}

		//Put it in the new vector
		if (create) {
			if (c_operator(&current_state, index)) can_projected = true;
		}
		else {
			if (c_dag_operator(&current_state, index)) can_projected = true;
		}

		if (can_projected) {
			sType index_vec;
			if (projected_states->where_is_element(current_state, &index_vec))
			{
				projected_vector[index_vec] = (T)phase * initial_vector[i];
			}
		}
	}
}

template <class T> void q_vector_creation(T* q_vector, T* prod_c_omega, 
                                int sspace_vec_size, T* subspace_evectors) {
	if(verbose < -9) std::cout << "Q_vector_creation(...)";
	/***************************************************************
	* Calculates the Q row (of the Q matrix) corresponding to the c_mu
	* 
	* Parameters
	* ----------
	* q_vector			: (T*) Array of the Q values linked to c_mu
	* prod_c_omega		: (T*) Fundamental state vector applied on c_mu 
    *                          (calculated before entering the function)
	* sspace_vec_size   : (int) Size of the subspace 
	* subspace_evectors : (T*) the eigen vectors of the subSpace matrix 
	*
	* Templates:
	* ----------
	* T : double, std::complex<double>
	*
	* Returns
	* -------
	* NONE
	*****************************************************************/
	//Compute <Omega|c_mu|phi_i>

	for (int n = 0; n < sspace_vec_size; n++){
		q_vector[n] = 0;
		for (int i = 0; i < sspace_vec_size; i++){
			q_vector[n] += subspace_evectors[i + n * sspace_vec_size]
                            * prod_c_omega[i];
		}
	}
}

template <class T, class StatesArrType> std::vector<double> compute_q_matrix (
        std::vector<T>* q_matrix, bool creation, int spin, 
        T* fund_state, StatesArrType* states_array) {
	/***************************************************************
	* Calculates the Q matrix
	* 
	* Parameters
	* ----------
	* q_matrix		: (std::vector<T>*) Q-Matrix
	* creation		: (bool) create or destroy green qM
	* spin			: (int) spin to create/destroy (1=up, 0=down)
	* fund_state		: (T*) fundamental vector |Omega>
	* states_array	: (StatesArrType) array of the states in the subspace
	*
	* Templates:
	* ----------
	* T             : double, std::complex<double>
    * StatesArrType : StatesR_T, StatesR_H
	*
	* Returns
	* -------
	* NONE
	*****************************************************************/
	if (verbose < -9) std::cout << "computeQ_matrix(...)" << std::endl;

	int sites = states_array->sys_hubP.n_sites;
	///Declare variables before parallization
	unsigned long new_space_len;
	uInt BL_space_size = BAND_LANCZOS_MAX_ITERATIONS;
	T* arr_BL; 
	std::vector<T> vec_BL;
	std::vector<T> prod_c_omega;
	//Eigen values of the BL subspace
	std::vector<double> BL_space_evalues;

	//Eigen vectors of the BL subspace
	std::vector<T> BL_space_evectors;

	//Projected excited states 
	StatesArrType* states_excited = states_array->clone();
	//Modify the number of  electrons
	states_excited->electrons.up += ((int)creation * 2 - 1) * spin ;
	states_excited->electrons.down += ((int)creation * 2 - 1) * (1-spin);

	//Projected Space
	green_space_projection(states_array, spin, creation, states_excited);
	
	//Number of times H is applied to generate new states
	if (verbose > 9) std::cout << "Before H excitation : " 
                               << states_excited->get_length() << std::endl;

	states_excited->subspace_condition_expanding();
	
	if (verbose > 9) std::cout << "After H excitation : " 
                               << states_excited->get_length() << std::endl;

	new_space_len = states_excited->get_length();
	arr_BL = new T[new_space_len * sites]();

	//Countains all the initial vectors for the band Lanczos algorithm
	//Creation of the p vectors for bandLanczos
	for (int i = 0; i < sites; i++){
		excited_vector_projection(creation, i, spin, fund_state, states_array,
                                  states_excited, arr_BL + i*new_space_len);
	}


	vec_BL = std::vector<T>(arr_BL, arr_BL + new_space_len * sites);
	//Band lanczos
	LanczosSolver<T,StatesArrType> LS;	
	BL_space_evalues = LS.band_lanczos_algorithm(
        &vec_BL, sites, new_space_len, states_excited, &BL_space_size, 
        &BL_space_evectors, &prod_c_omega);

	//Clear mem
	std::vector<T>().swap(vec_BL);
	states_excited->remove_all();
	delete[] arr_BL;
	
	//Repeat for each nu/mu pair
	//All Green values vector
	if(verbose > 9) std::cout << "Q MATRIX STARTING" << std::endl;
	*q_matrix = std::vector<T>(sites * BL_space_size, 0);

	for (int j = 0; j < sites; j++) {
		q_vector_creation(q_matrix->data() + BL_space_size * j, 
                    prod_c_omega.data() + j * BL_space_size, BL_space_size,
                    BL_space_evectors.data());

	}
	delete states_excited;

	return BL_space_evalues;
}

//Greens
template <class T, class StatesArrType> void compute_q_matrix_band_lanczos(
        int spin, T* fund_state, double fundE, 
        StatesArrType* const states_array, greenParam gP,int deg) { 
	/***************************************************************
	* Calculates the Q matrix and writes them in a qMatrices.txt
	* 
	* Parameters
	* ----------
	* spin			: (int) spin to create/destroy (1=up, 0=down)
	* fund_state	: (T*) fundamental vector |Omega>
	* fundE			: (double) fundamental energy
	* states_array	: (StatesArrType*) array of the states in the subspace
	* gP			: (greenParam) green params
	* deg			: (int) degeneracy of the fundamental
	*
	* Templates:
	* ----------
	* T : double, std::complex<double>
    * StatesArrType : StatesR_T, StatesR_H
	*
	* Returns
	* -------
	* NONE
	*****************************************************************/
	if (verbose < -4) std::cout << "compute_q_matrix_band_lanczos(...)" 
                                << std::endl;

	int size = states_array->get_length();
	int sites = states_array->sys_hubP.n_sites;

	//////Wrtiting the Q-Matrix and the eigen values in a txt file
	Electrons elec = states_array->electrons;
	float perc_used = (float)states_array->get_length() / 
                        (comb(sites, elec.up) * comb(sites, elec.down));
	std::string written_q_matrix = write_q_system(fundE, (double)gP.g_eta, elec,
                                            &states_array->sys_hubP, perc_used);

	written_q_matrix += "\n# Eigen values E -- Q-Matrixes E\n";

	for (int d = 0; d < deg; d++) {
		std::vector<T> q_matrix_e;
		std::vector<double> BL_space_evalues_e = compute_q_matrix(
            &q_matrix_e, true, spin, fund_state + d*size, states_array);
		written_q_matrix += write_q_matrix_string(
            &q_matrix_e, &BL_space_evalues_e, sites)+"\n";
	}

	written_q_matrix += "\n# Eigen values H -- Q-Matrixes H\n";
	for (int d = 0; d < deg; d++) {
		std::vector<T> q_matrix_h;
		std::vector<double> BL_space_evalues_h = compute_q_matrix(
            &q_matrix_h, false, spin, fund_state + d*size, states_array);
		written_q_matrix += write_q_matrix_string(
            &q_matrix_h, &BL_space_evalues_h, sites)+"\n";
	}
	
	char cwd[PATH_MAX];
	if(!getcwd(cwd, sizeof(cwd))) std::cout << "Problem occured getting CWD"
                                            << std::endl;
	std::string txt_name = "/qMatrices.txt";
	const std::string out_file_name = cwd + txt_name;
	std::ofstream out_file;
	//Write Q-matrix and ev
	out_file.open(out_file_name);
	out_file << written_q_matrix;
	out_file.close();
	//////
}

template <class StatesArrType> void compute_green_long(
        int spin, std::vector<double>* fund_state, double fundE, 
        StatesArrType* const states_array, greenParam gP, int deg){
	/***************************************************************
	* Calculates the Q matrix and writes them in a qMatrices.txt
	* 
	* Parameters
	* ----------
	* spin			: (int) spin to create/destroy (1=up, 0=down)
	* fund_state	: (std::vector<T>*) fundamental vector |Omega>
	* fundE			: (double) fundamental energy
	* states_array	: (StatesArrType*) array of the states in the subspace
	* gP			: (greenParam) green params
	* deg			: (int) degeneracy of the fundamental
	*
	* Templates:
	* ----------
	* StatesArrType : StatesR_T, StatesR_H
	*
	* Returns
	* -------
	* NONE
	*****************************************************************/
	if(verbose < -4) std::cout<<"computeGreen_long(...)"<<std::endl;

	int sites = states_array->sys_hubP.n_sites;
	//Projected excited states 
	StatesArrType* states_excited_e = states_array->clone();
	states_excited_e->electrons.up += spin ;
	states_excited_e->electrons.down += (1-spin);
	StatesArrType* states_excited_h = states_array->clone();
	states_excited_h->electrons.up -= spin ;
	states_excited_h->electrons.down -= (1-spin);
	
	//Projected Space
	green_space_projection(states_array, spin, true, states_excited_e);
	green_space_projection(states_array, spin, false, states_excited_h);

	states_excited_e->subspace_condition_expanding();
	states_excited_h->subspace_condition_expanding();

	if (verbose > 99) {
		std::cout<<"PROJ STATES E"<<std::endl;
		states_excited_e->show_all_states();
		std::cout<<"PROJ STATES H"<<std::endl;
		states_excited_h->show_all_states();
	}

	//Wrtiting the Q-Matrix and the eigen values in a txt file
	Electrons elec = states_array->electrons;
	float perc_used = (float)states_array->get_length() / 
                            (comb(sites, elec.up) * comb(sites, elec.down));
	std::string written_q_matrix = write_q_system(fundE, (double)gP.g_eta, elec,
                                            &states_array->sys_hubP, perc_used);
	
	//Vectors
	//ELECTONS
	int new_space_len_e = states_excited_e->get_length();

	written_q_matrix += "\n# Eigen values E -- Q-Matrixes E\n";

	for (int m = 0; m < deg; m++) {
		std::vector<double> q_matrix_e = std::vector<double>(
                                                    sites*new_space_len_e, 0);
		std::vector<double> eigen_e;
		if (new_space_len_e > 0) {
			double* arr_BL_e = new double[new_space_len_e * sites]();
			//Creation of the vectors c_mu^(dag)|Omega>
			for (int i = 0; i < sites; i++){
				excited_vector_projection(
                    true, i, spin, 
                    fund_state->data() + states_array->get_length() * m,
                    states_array, states_excited_e, 
                    arr_BL_e + i * new_space_len_e);
			}

			//Hamiltonian matrices
			double* hE = new double[new_space_len_e*new_space_len_e]();
			states_excited_e->matrix_creation(hE);

			char jobs = 'V', uplo='U';
			double* eigen_value_e = new double[new_space_len_e]();
			int lwork = new_space_len_e*(new_space_len_e+1);
			double* work = new double[lwork];
			int info;

			//Eigen values of hE|E> = E|E>
			dsyev_(&jobs, &uplo, &new_space_len_e, hE, &new_space_len_e,
                   eigen_value_e, work, &lwork, &info);

			delete[] work; 

			eigen_e = std::vector<double>(eigen_value_e, 
                                            eigen_value_e + new_space_len_e);
			//<OMEGA|c Ue
            char trans_a = 'N', trans_b = 'T';
			dgemm_(&trans_a, &trans_b, &sites, &new_space_len_e, &new_space_len_e, &ALPHA_D, 
                  arr_BL_e, &new_space_len_e, hE, &new_space_len_e, &BETA_D,
                  q_matrix_e.data(), &new_space_len_e);

			delete[] hE; delete[] arr_BL_e;	delete[] eigen_value_e; 	
		}
		written_q_matrix += write_q_matrix_string(
                                            &q_matrix_e, &eigen_e, sites)+"\n";
	}
	delete states_excited_e;


	//HOLES
	int new_space_len_h = states_excited_h->get_length();


	//Wrtiting the Q-Matrix and the eigen values in a txt file
	written_q_matrix += "\n# Eigen values H -- Q-Matrixes H\n";

	for (int m = 0; m < deg; m++) {
		std::vector<double> q_matrix_h = std::vector<double>(
                                                    sites*new_space_len_h, 0);
		std::vector<double> eigen_h;
		if (new_space_len_h > 0) {
			double* arr_BL_h = new double[new_space_len_h * sites]();
			//Creation of the vectors c_mu^(dag)|Omega>
			for (int i = 0; i < sites; i++){
				excited_vector_projection(
                    false, i, spin, 
                    fund_state->data() + states_array->get_length() * m, 
                    states_array, states_excited_h, 
                    arr_BL_h + i * new_space_len_h);
			}

			//Hamiltonian matrices
			double* hH = new double[new_space_len_h*new_space_len_h]();
			states_excited_h->matrix_creation(hH);

			char jobs = 'V', uplo='U';
			double* eigen_value_h = new double[new_space_len_h]();
			int lwork = new_space_len_h*(new_space_len_h+1);
			double* work = new double[lwork];
			int info;

			//Eigen values of hH|E> = E |E>
			dsyev_(&jobs, &uplo, &new_space_len_h, hH, &new_space_len_h, 
                   eigen_value_h, work, &lwork, &info);

			delete[] work; 
			
			eigen_h = std::vector<double>(eigen_value_h, 
                                            eigen_value_h + new_space_len_h);
			//<OMEGA|c Uh
            char trans_a = 'N', trans_b = 'T';
			dgemm_(&trans_a, &trans_b, &sites, &new_space_len_h, &new_space_len_h, &ALPHA_D, 
                  arr_BL_h, &new_space_len_h, hH, &new_space_len_h, &BETA_D, 
                  q_matrix_h.data(), &new_space_len_h);

			delete[] hH; delete[] arr_BL_h; delete[] eigen_value_h; 
		}
		written_q_matrix += write_q_matrix_string(&q_matrix_h, 
                                                    &eigen_h, sites)+"\n";
	}
	delete states_excited_h;

	char cwd[PATH_MAX];
	if(!getcwd(cwd, sizeof(cwd))) std::cout << "Problem occured getting CWD"
                                            << std::endl;
	std::string txt_name = "/qMatrices.txt";
	const std::string out_file_name = cwd + txt_name;
	std::ofstream out_file;
	//Write Q-matrix and ev
	out_file.open(out_file_name);
	out_file << written_q_matrix;
	out_file.close();
	//////
	
}
#endif

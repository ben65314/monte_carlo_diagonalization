#include "LanczosSolver.h"

#ifndef __greenFunctions_h__
#define __greenFunctions_h__

template <class StatesArrType> void c_subSpace(StatesArrType* const sub_space, int location, int spin, bool creation, StatesArrType* c_sub_space){
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
	unsigned long len = sub_space->getLength();
	
	int index = sites - location - 1 + spin * sites;
	if (creation) {
		for (unsigned long i = 0; i < len; i++) {
			sType temp = sub_space->getAt(i);
			if (c_operator(&temp, index)) { 
				c_sub_space->add(temp);
			}
		}
	}
	else {
		for (unsigned long i = 0; i < len; i++) {
			sType temp = sub_space->getAt(i);
			if (c_dag_operator(&temp, index)) { 
				c_sub_space->add(temp);
			}
		}
	}
}

template <class StatesArrType> void greenSpaceProjection(StatesArrType* const origin_sub_space, int spin, bool creation, StatesArrType* proj_sub_space){
	/****************************************************************
	* Applies all the possible c_mu to create all the states of the projected space
	* 
	* Parameter
	* ----------
	* origin_sub_space	: (StatesArrType) States who will be applied the c operators
	* spin				: (int) value of the spin created or destroyed (1=up, 0=down)
	* creation			: (bool) Is c a creation or annihilation operator
	* proj_sub_space	: (StatesArrType) States resulting of the c appliacation
	*
	* Returns
	* -------
	* NONE
	*****************************************************************/
	if (verbose < -9) std::cout<<"greenSpaceProjection(...)"<<std::endl;
	//originSubSpace->showAllStates();
	unsigned int sites = origin_sub_space->sys_hubP.n_sites;
	proj_sub_space->removeAll();
	//StatesArr* temp;
	for (unsigned int i = 0; i < sites; i++) {
		c_subSpace(origin_sub_space, i, spin, creation, proj_sub_space);
		//projSubSpace->merge(temp);
	}
}

std::string writeQ_System(double fundE, double eta, Electrons elec, hubbardParam* hubP_to_write, float percStates){
	/****************************************************************
	* Writes the system parameters for the qMatrices.txt file
	* 
	* Parameter
	* ----------
	* fundE			 : (double) Fundamental Energy
	* eta			 : (double) Width of the Lorentz
	* elec			 : (Electrons) Electrons of the bloc system 
	* hubP_to_write: (hubbardParam*) System parameters to write
	* perc_states	 : (float) {ercentage of states taken
	*
	* Templates
	* ---------
	* T				 : float, double
	*
	* Returns
	* -------
	* text			 : (std::string) String of parameters
	*****************************************************************/
	std::string text = "";
	text += to_string_p(hubP_to_write->n_sites,0) + "\t";
	text += to_string_p(hubP_to_write->u,0) + "\t";
	text += to_string_p(hubP_to_write->mu,0) + "\t";
	text += to_string_p((elec.up+elec.down),0) + "\t";
	text += to_string_p((elec.up-elec.down),0) + "\t";
	text += to_string_p(fundE,16) + "\t";
	text += to_string_p(eta,5) + "\t";
	text += to_string_p(percStates,2) + "\t";
	return text;
}

template <class T> std::string writeQ_MatrixString(std::vector<T>* QM, std::vector<double>* eigen_val, int ns) {
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
			text += to_string_pq(QM->at(i * eigen_val->size() + j),4,14) + "\t";
		}
		if (j < eigen_val->size()-1) text+="\n";
	}

	return text;
}

template <class T> std::string writeCF_MatrixString(std::vector<double>* all_alpha, std::vector<T>* all_beta, std::vector<int>* len_param) {
	int max = *std::max_element(len_param->begin(), len_param->end());

	std::string text = "!a\n";
	for (uInt j = 0; j < max; j++) {
		int iterate = 0 ;
		for (uInt i = 0; i < len_param->size(); i++) {
			if (j < len_param->at(i)) {
				text += to_string_pq(all_alpha->at(j+iterate),4,14) + "\t";
			}
			else {
				text += to_string_pq(0.0,4,14) + "\t";
				
			}
			iterate += len_param->at(i);
		}
		text+="\n";
	}
	text+="!b\n";

	for (uInt j = 0; j < max+1; j++) {
		int iterate = 0;
		for (uInt i = 0; i < len_param->size(); i++) {
			if (j < len_param->at(i)) {
				text += to_string_pq(all_beta->at(j+iterate),4,14) + "\t";
			}
			else {
				text += to_string_pq(0.0,4,14) + "\t";
			}
			iterate += len_param->at(i)+1;
		}
		text+="\n";
	}

	return text;
}

template <class T, class StatesArrType> void excitedVectorProjection(bool create, unsigned int site, int spin, const T* initial_vector, const StatesArrType* initial_states, StatesArrType* projected_states, T* projected_vector) {
	/****************************************************************
	* Transform a vector countained in a space to the according vector in the excited space. 
	* 
	* Parameters
	* ----------
	* create			: (bool) Are we creating or annihilating an electron to get to the excited space
	* site				: (int) which site is the electron or hole created
	* spin				: (int) value of the spin created or destroyed (1=up, 0=down)
	* initial_vector	: (const T*) the initial vector in the initialStates space
	* initial_states	: (const StatesArr*) the initial state space
	* projected_states	: (StatesArr*) the projected state space
	* projected_vector	: (T*) the projected vector in the projectedStates space
	* 
	* Templates:
	* ----------
	* T					: double, std::complex<double>
	*
	* Returns
	* -------
	* NONE
	*****************************************************************/
	sType current_state;
	//printVector(initialVector,initialStates->getLength());
	int ns = initial_states->sys_hubP.n_sites; //number of sites
	int index = ns - site - 1 + spin * ns;
	
	for(long unsigned i = 0; i<initial_states->getLength();i++){
		//Exciting the state
		current_state = initial_states->getAt(i);
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
			if (projected_states->whereIsElement(current_state,&index_vec))
			{
				projected_vector[index_vec] = (T)phase * initial_vector[i];
			}
		}
	}
}

template <class T> void Q_vector_creation(T* q_vector, T* prod_c_omega, int sspace_vec_size, T* subspace_Evectors) {
	if(verbose < -9) std::cout<<"Q_vector_creation(...)";
	/***************************************************************
	* Calculates the Q row (of the Q matrix) corresponding to the c_mu
	* 
	* Parameters
	* ----------
	* q_vector			: (T*) Array of the Q values linked to c_mu
	* prod_c_omega		: (T*) Fundamental state vector applied on c_mu (calculated before entering the function)
	* sspace_vec_size		: (int) Size of the subspace 
	* subspace_Evectors : (T*) the eigen vectors of the subSpace matrix 
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
			q_vector[n] += subspace_Evectors[i + n * sspace_vec_size] * prod_c_omega[i];
		}
	}
}

template <class T, class StatesArrType> std::vector<double> computeQ_matrix (std::vector<T>* Q_matrix, bool creation, int spin, T* fundState, StatesArrType* states_array) {
	/***************************************************************
	* Calculates the Q matrix
	* 
	* Parameters
	* ----------
	* Q_matrix		: (std::vector<T>*) Q-Matrix
	* creation		: (bool) create or destroy green qM
	* spin			: (int) spin to create/destroy (1=up, 0=down)
	* fundState		: (T*) fundamental vector |Omega>
	* statesArray	: (StatesArrType) array of the states in the subspace
	*
	* Templates:
	* ----------
	* T : double, std::complex<double>
	*
	* Returns
	* -------
	* NONE
	*****************************************************************/
	if(verbose < -9) std::cout<<"computeQ_matrix(...)"<<std::endl;

	int sites = states_array->sys_hubP.n_sites;
	///Declare variables before parallization
	unsigned long newSpaceLen;
	uInt BLspaceSize = BAND_LANCZOS_MAX_ITERATIONS;
	T* arr_BL; 
	std::vector<T> vecBL;
	std::vector<T> ProdCOmega;
	//Eigen values of the BL subspace
	std::vector<double> BLspace_Evalues;

	//Eigen vectors of the BL subspace
	std::vector<T> BLspace_Evectors;

	//Projected excited states 
	StatesArrType* statesExcited = states_array->clone();
	//Modify the number of  electrons
	statesExcited->electrons.up += ((int)creation * 2 - 1) * spin ;
	statesExcited->electrons.down += ((int)creation * 2 - 1) * (1-spin);

	//Projected Space
	greenSpaceProjection(states_array, spin, creation, statesExcited);
	
	//Number of times H is applied to generate new states
	if (verbose > 9) std::cout<<"Before H excitation : "<<statesExcited->getLength()<<std::endl;

	statesExcited->subspace_condition_expanding();
	
	if (verbose > 9) std::cout<<"After H excitation : "<<statesExcited->getLength()<<std::endl;

	newSpaceLen = statesExcited->getLength();
	arr_BL = new T[newSpaceLen * sites]();

	//Countains all the initial vectors for the band Lanczos algorithm
	//Creation of the p vectors for bandLanczos
	for (int i = 0; i < sites; i++){
		excitedVectorProjection(creation, i, spin, fundState, states_array, statesExcited, arr_BL + i * newSpaceLen);
	}


	vecBL = std::vector<T>(arr_BL, arr_BL + newSpaceLen * sites);
	//Band lanczos
	LanczosSolver<T,StatesArrType> LS;	
	BLspace_Evalues = LS.band_lanczos_algorithm(&vecBL, sites, newSpaceLen, statesExcited, &BLspaceSize, &BLspace_Evectors, &ProdCOmega);

	//printMatrix(BLspace_Evectors.data(),)
	//Clear mem
	std::vector<T>().swap(vecBL);
	statesExcited->removeAll();
	delete[] arr_BL;
	
	//Repeat for each nu/mu pair
	//All Green values vector
	if(verbose > 9) std::cout<<"Q MATRIX STARTING"<<std::endl;
	*Q_matrix = std::vector<T>(sites * BLspaceSize,0);

	for (int j = 0; j < sites; j++) {
		Q_vector_creation(Q_matrix->data() + BLspaceSize * j, ProdCOmega.data() + j * BLspaceSize, BLspaceSize, BLspace_Evectors.data());

	}
	delete statesExcited;

	return BLspace_Evalues;
}

//Greens
template <class T, class StatesArrType> void computeQMatrixBandLanczos(int spin, T* fundState, double fundE, StatesArrType* const states_array, greenParam gP,int deg) { 
	/***************************************************************
	* Calculates the Q matrix and writes them in a qMatrices.txt
	* 
	* Parameters
	* ----------
	* spin			: (int) spin to create/destroy (1=up, 0=down)
	* fundState		: (T*) fundamental vector |Omega>
	* fundE			: (double) fundamental energy
	* statesArray	: (StatesArrType*) array of the states in the subspace
	* gP			: (greenParam) green params
	* deg			: (int) degeneracy of the fundamental
	*
	* Templates:
	* ----------
	* T : double, std::complex<double>
	*
	* Returns
	* -------
	* NONE
	*****************************************************************/
	if(verbose < -4) std::cout<<"computeQMatrixBandLanczos(...)"<<std::endl;

	int size = states_array->getLength();
	int sites = states_array->sys_hubP.n_sites;

	//////Wrtiting the Q-Matrix and the eigen values in a txt file
	Electrons elec = states_array->electrons;
	float percUsed = (float)states_array->getLength() / (comb(sites,elec.up) * comb(sites,elec.down));
	std::string writtenQMatrix = writeQ_System(fundE, (double)gP.g_eta, elec, &states_array->sys_hubP, percUsed);

	writtenQMatrix += "\n# Eigen values E -- Q-Matrixes E\n";

	for (int d = 0; d < deg; d++) {
		std::vector<T> Q_matrix_e;
		std::vector<double> BLspace_Evalues_e = computeQ_matrix (&Q_matrix_e, true, spin, fundState + d*size, states_array);
		writtenQMatrix += writeQ_MatrixString(&Q_matrix_e, &BLspace_Evalues_e, sites)+"\n";
	}

	writtenQMatrix += "\n# Eigen values H -- Q-Matrixes H\n";
	for (int d = 0; d < deg; d++) {
		std::vector<T> Q_matrix_h;
		std::vector<double> BLspace_Evalues_h = computeQ_matrix (&Q_matrix_h, false, spin, fundState + d*size, states_array);
		writtenQMatrix += writeQ_MatrixString(&Q_matrix_h, &BLspace_Evalues_h, sites)+"\n";
	}
	
	char cwd[PATH_MAX];
	if(!getcwd(cwd, sizeof(cwd))) std::cout<<"Problem occured getting CWD"<<std::endl;
	std::string txtName = "/qMatrices.txt";
	const std::string outFileName = cwd + txtName;
	std::ofstream outFile;
	//Write Q-matrix and ev
	outFile.open(outFileName);
	outFile<<writtenQMatrix;
	outFile.close();
	//////
}

template <class StatesArrType> void computeGreen_long(int spin, std::vector<double>* fundState, double fundE, StatesArrType* const states_array, greenParam gP, int deg){
	/***************************************************************
	* Calculates the Q matrix and writes them in a qMatrices.txt
	* 
	* Parameters
	* ----------
	* spin			: (int) spin to create/destroy (1=up, 0=down)
	* fundState		: (std::vector<T>*) fundamental vector |Omega>
	* fundE			: (double) fundamental energy
	* states_array	: (StatesArrType*) array of the states in the subspace
	* gP			: (greenParam) green params
	* deg			: (int) degeneracy of the fundamental
	*
	* Templates:
	* ----------
	* T : double, std::complex<double>
	*
	* Returns
	* -------
	* NONE
	*****************************************************************/
	if(verbose < -4) std::cout<<"computeGreen_long(...)"<<std::endl;

	int sites = states_array->sys_hubP.n_sites;
	//Projected excited states 
	StatesArrType* statesExcited_E = states_array->clone();
	statesExcited_E->electrons.up += spin ;
	statesExcited_E->electrons.down += (1-spin);
	StatesArrType* statesExcited_H = states_array->clone();
	statesExcited_H->electrons.up -= spin ;
	statesExcited_H->electrons.down -= (1-spin);
	
	//Projected Space
	greenSpaceProjection(states_array, spin, true, statesExcited_E);
	greenSpaceProjection(states_array, spin, false, statesExcited_H);

	statesExcited_E->subspace_condition_expanding();
	statesExcited_H->subspace_condition_expanding();

	if (verbose > 99) {
		std::cout<<"PROJ STATES E"<<std::endl;
		statesExcited_E->showAllStates();
		std::cout<<"PROJ STATES H"<<std::endl;
		statesExcited_H->showAllStates();
	}

	//////Wrtiting the Q-Matrix and the eigen values in a txt file
	Electrons elec = states_array->electrons;
	float percUsed = (float)states_array->getLength() / (comb(sites,elec.up) * comb(sites,elec.down));
	std::string writtenQMatrix = writeQ_System(fundE, (double)gP.g_eta, elec, &states_array->sys_hubP, percUsed);
	
	//Vectors
	//ELECTONS
	int newSpaceLen_E = statesExcited_E->getLength();

	writtenQMatrix += "\n# Eigen values E -- Q-Matrixes E\n";

	for (int m = 0; m < deg; m++) {
		std::vector<double> Q_matrix_E = std::vector<double>(sites*newSpaceLen_E,0);
		std::vector<double> eigen_E;
		if (newSpaceLen_E > 0) {
			double* arr_BL_E = new double[newSpaceLen_E * sites]();
			//Creation of the vectors c_mu^(dag)|Omega>
			for (int i = 0; i < sites; i++){
				excitedVectorProjection(true, i, spin, fundState->data() + states_array->getLength() * m, states_array, statesExcited_E, arr_BL_E + i * newSpaceLen_E);
			}

			//Hamiltonian matrices
			double* hE = new double[newSpaceLen_E*newSpaceLen_E]();
			statesExcited_E->matrixCreation(hE);

			char jobs = 'V', uplo='U';
			double* eigen_EE = new double[newSpaceLen_E]();
			int lwork = newSpaceLen_E*(newSpaceLen_E+1);
			double* work = new double[lwork];
			int info;

			//Eigen values of hE|E> = E|E>
			dsyev_(&jobs,&uplo,&newSpaceLen_E,hE,&newSpaceLen_E,eigen_EE,work,&lwork,&info,1,1);

			delete[] work; 

			eigen_E = std::vector<double>(eigen_EE, eigen_EE + newSpaceLen_E);
			//<OMEGA|c Ue
			cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,sites,newSpaceLen_E,newSpaceLen_E,ALPHA_D,arr_BL_E,newSpaceLen_E,hE,newSpaceLen_E,BETA_D,Q_matrix_E.data(),newSpaceLen_E);

			delete[] hE; delete[] arr_BL_E;	delete[] eigen_EE; 	
		}
		writtenQMatrix += writeQ_MatrixString(&Q_matrix_E, &eigen_E, sites)+"\n";
	}
	delete statesExcited_E;


	//HOLES
	int newSpaceLen_H = statesExcited_H->getLength();


	//////Wrtiting the Q-Matrix and the eigen values in a txt file
	writtenQMatrix += "\n# Eigen values H -- Q-Matrixes H\n";

	for (int m = 0; m < deg; m++) {
		std::vector<double> Q_matrix_H = std::vector<double>(sites*newSpaceLen_H,0);
		std::vector<double> eigen_H;
		if (newSpaceLen_H > 0) {
			double* arr_BL_H = new double[newSpaceLen_H * sites]();
			//Creation of the vectors c_mu^(dag)|Omega>
			for (int i = 0; i < sites; i++){
				excitedVectorProjection(false, i, spin, fundState->data() + states_array->getLength() * m, states_array, statesExcited_H, arr_BL_H + i * newSpaceLen_H);
			}

			//Hamiltonian matrices
			double* hH = new double[newSpaceLen_H*newSpaceLen_H]();
			statesExcited_H->matrixCreation(hH);

			//std::cout<<"MATRIX HOLES:"<<std::endl;
			//printMatrix(hH,newSpaceLen_H,newSpaceLen_H,3,0);

			char jobs = 'V', uplo='U';
			double* eigen_HH = new double[newSpaceLen_H]();
			int lwork = newSpaceLen_H*(newSpaceLen_H+1);
			double* work = new double[lwork];
			int info;

			//Eigen values of hH|E> = E |E>
			dsyev_(&jobs,&uplo,&newSpaceLen_H,hH,&newSpaceLen_H,eigen_HH,work,&lwork,&info,1,1);

			delete[] work; 
			
			//std::cout<<"H|OMEGA>:"<<std::endl;
			//printMatrix(arr_BL_H,(int)sites,newSpaceLen_H,3,6);
			//std::cout<<"H vectors"<<std::endl;
			//printMatrix(hH,newSpaceLen_H,newSpaceLen_H,5,5);

			eigen_H = std::vector<double>(eigen_HH, eigen_HH + newSpaceLen_H);
			//<OMEGA|c Uh
			cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,sites,newSpaceLen_H,newSpaceLen_H,ALPHA_D,arr_BL_H,newSpaceLen_H,hH,newSpaceLen_H,BETA_D,Q_matrix_H.data(),newSpaceLen_H);

			delete[] hH; delete[] arr_BL_H; delete[] eigen_HH; 
		}
		writtenQMatrix += writeQ_MatrixString(&Q_matrix_H, &eigen_H, sites)+"\n";
	}
	delete statesExcited_H;

	char cwd[PATH_MAX];
	if(!getcwd(cwd, sizeof(cwd))) std::cout<<"Problem occured getting CWD"<<std::endl;
	std::string txtName = "/qMatrices.txt";
	const std::string outFileName = cwd + txtName;
	std::ofstream outFile;
	//Write Q-matrix and ev
	outFile.open(outFileName);
	outFile<<writtenQMatrix;
	outFile.close();
	//////
	
}
#endif

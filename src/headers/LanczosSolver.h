#pragma once

#include "ModulesStates/StatesR_T.h"
#include "ModulesStates/StatesR_H.h"


//GENERIC TEMPLATE
	template<class T, class StatesArrType> class LanczosSolver;

//DOUBLE TEMPLATE
template<class StatesArrType> class LanczosSolver<double,StatesArrType>{
	public:
	//The three functions are used to compute the fund energy while keeping only three vectors in memory
	void LanczosEnergy(std::vector<double>* fundState_lanczos_basis, double* init_vector, StatesArrType* sArr, std::vector<double>* alpha, std::vector<double>* beta, double* fund_energy, int* iter, int* deg, double epsilon){
		/*******************************************************
		* Computes the alpha/beta of the tridiagonal matrix and find the converged fund energy and fund vector in the reduced space.
		*
		* Parameters
		* ----------
		* fund_state_lanczos_basis	: (std::vector<double>*) fundamental eigen vector in the Lanczos basis
		* init_vector				: (double) initial vector |phi_0>
		* sArr						: (StatesArrType) array of states object
		* alpha						: (std::vector<double>*) alphas
		* beta						: (std::vector<double>*) betas
		* fund_energy				: (double*) fundamental energy receptacle
		* iter						: (int*) number of iterations done
		* deg						: (int*) degeneracy of the fund vector
		* epsilon					: (double) convergence precision
		*
		* Returns
		* -------
		* NONE
		*******************************************************/
		//Sets size of alpha and beta
		int size = sArr->getLength();
		alpha->clear();	alpha->reserve(*iter);
		beta->clear();	beta->reserve(*iter);

		//Two main vectors
		std::vector<double> r(init_vector, init_vector + size);
		std::vector<double> q(size);

		//Energies to converge
		double prev_iter_energy = 1000;
		double energy = 100;
		
		int current_iteration = 1;

		bool converged = false;

		while (!converged) {
			if(beta->size()) {
				cblas_dscal(size, 1/beta->back(), r.data(), 1);
				cblas_dscal(size, -beta->back(), q.data(), 1);
			}

			std::vector<double> H_tmp(size);
			//Applies the vector r on the matrix H and stores it in H_tmp
			//std::cout<<"r"<<std::endl;
			//printVector(r.data(),r.size());
			sArr->H(H_tmp.data(),r.data());
			//std::cout<<"Hr"<<std::endl;
			//printVector(H_tmp.data(),H_tmp.size());

			cblas_daxpy(size,1,H_tmp.data(),1,q.data(),1);	//q = q + H*r
			//std::cout<<"q"<<std::endl;
			//printVector(q.data(),q.size());
			//swap q <-> r
			cblas_dswap(size,q.data(),1,r.data(),1);

			double dot_product = cblas_ddot(size,q.data(),1,r.data(),1);
			//std::cout<<"DOT PROD = "<<dotProd<<std::endl;
			alpha->push_back(dot_product);

			cblas_daxpy(size,-alpha->back(),q.data(),1,r.data(),1);	//r = r - q*alpha

			beta->push_back(cblas_dnrm2(size,r.data(),1));

			//Arrays for tridiag solve
			double* arr_a = new double[alpha->size()];
			double* arr_b = new double[beta->size()];
			std::copy(alpha->begin(),alpha->end(),arr_a);
			std::copy(beta->begin(),beta->end(),arr_b);
			
			//Parameters for solver
			char jobs = 'V'; 
			int n = current_iteration; 
			int info;
			double tt = 1.0;
			double* vecs = new double[n*n];
			double* work = new double[2*n];

			//Solving Energy
			dstev_(&jobs,&n,arr_a,arr_b,vecs,&n,work,&info,tt);
			//Fund energy

			prev_iter_energy = energy;
			energy = arr_a[0];


			if (abs(prev_iter_energy - energy) < epsilon && current_iteration > 3) {
				converged = true;
				*deg = degFundamentalCheck(arr_a,n);
				fundState_lanczos_basis->clear();
				*fundState_lanczos_basis = std::vector<double>(vecs, vecs + n*(*deg));
			}

			delete[] vecs; delete[] work; delete[]arr_a; delete[] arr_b;

			if(current_iteration%10==0 && verbose > 4){
				std::cout << "Lanczos Energy current iteration :" << current_iteration << ":"<<abs(prev_iter_energy - energy)<< std::endl;
			}
			if(current_iteration == size) {
				break;
			}
			current_iteration++;
			
		}
		if(verbose > 4){
			std::cout << "Lanczos iteration used:" << current_iteration << std::endl;
		}
		*iter = current_iteration;
		*fund_energy = energy;
	}
	void LanczosVectors(std::vector<double>* fundState_lanczos_basis, double* gs, StatesArrType* sArr, std::vector<double>* alpha, std::vector<double>* beta, int* deg){
		/*******************************************************
		* Convert the fund vector in the reduced space to the original space.
		*
		* Parameters
		* ----------
		* fund_state_lanczos_basis	: (std::vector<double>*) fundamental eigen vector in the Lanczos basis
		* gs						: (double) initial vector |phi_0> and on out groundstate
		* sArr						: (StatesArrType) array of states object
		* alpha						: (std::vector<double>*) alphas
		* beta						: (std::vector<double>*) betas
		* deg						: (int*) degeneracy of the fund vector
		*
		* Returns
		* -------
		* NONE
		*******************************************************/
		int size = sArr->getLength();
		int size_proj = alpha->size();

		std::vector<double> r(gs, gs + size);
		std::vector<double> q(size);

		for (int d = 0; d < *deg; d++)
		cblas_dscal(size, fundState_lanczos_basis->at(0+size_proj*d), gs+d*size, 1);


		for (unsigned int j = 1; j <fundState_lanczos_basis->size(); j++) {
			std::vector<double> H_tmp(size);
			sArr->H(H_tmp.data(),r.data());

			cblas_daxpy(size,1,H_tmp.data(),1,q.data(),1);	//q = q + H*r

			cblas_daxpy(size,-alpha->at(j-1),r.data(),1,q.data(),1);	//r = r - q*alpha
			for (unsigned int i = 0; i < r.size(); i++) {
				double tmp = r[i];
				r[i] = q[i]/beta->at(j-1);
				q[i] = -beta->at(j-1)*tmp;
			}

			for (int d = 0; d < *deg; d++) {
				double scal = fundState_lanczos_basis->at(j+size_proj*d);
				cblas_daxpy(size,scal,r.data(),1,gs+d*size,1);	//r = r - q*alpha
			}
			if(j%10==0 && verbose > 4){
				std::cout << "Lanczos Vector current iteration :" << j << std::endl;
			}
		}
	}
	double lanczosAlgorithm(std::vector<double>* fundState, StatesArrType* sArr, int* deg, double epsilon = 10e-12) {
		/*************************************************
		* Redefines a given matrix with the Lanczos Algorithm without needing the Hamiltonian matrix
		*
		* Parameters
		* ----------
		* fundState	: (std::vector<double>*) fundamental state of the system
		* sArr		: (StatesArr*) States used in the subspace
		* deg			: (int*) degeneracy counter	 
		* epsilon		: (double) Convergence acceptability
		*
		* Returns
		* -------
		* currentEnergy : (double) fundamental energy of the system
		***************************************************/
		if(verbose == -1) std::cout<<"double lanczosAlgorithm(double...) called"<<std::endl;

		double fundEnergy;
		std::vector<double> alpha, beta;
		std::vector<double> fundState_lanczosBasis;
		int nIterations = 1000;

		// Random Initial Vector
		initialVector(sArr->getLength(),fundState->data());
		
		LanczosEnergy(&fundState_lanczosBasis,fundState->data(),sArr,&alpha,&beta,&fundEnergy,&nIterations,deg,epsilon);
		if (verbose > 9) std::cout<<"DEGENERACY:"<<*deg<<std::endl;
		//Increase the size of the fundState according to the degeneracy
		for (int i = 1; i < *deg; i++) {
			fundState->insert(fundState->end(),fundState->begin(),fundState->begin()+sArr->getLength());
		}
		LanczosVectors(&fundState_lanczosBasis,fundState->data(),sArr,&alpha,&beta,deg);

		return fundEnergy;
	}	


	std::vector<double> bandLanczosAlgorithm(std::vector<double>* vk, uInt n_bk, unsigned long len_bk, StatesArrType* sArr, uInt* nIter, std::vector<double>* sub_space_vectors, std::vector<double>* product_c_omega, double dtol = 10e-10){
		/*************************************************
		* Band lanczos algorithm
		*
		* Parameters
		* ----------
		* vk				: (std::vector<double>*) The n_bk initial vectors
		* n_bk				: (uInt) number of band vectors
		* len_bk			: (unsigned long) length of the vectors
		* sArr				: (StatesArrType*) states array
		* nIter				: (uInt*) band lanczos iteration done
		* sub_space_vectors	: (std::vector<double>* ) 
		* product_c_omega	: (std::vector<double>*) <Omega|c|phi>
		* dtol				: (double) convergence tolerance
		*
		*
		* Returns
		* -------
		* energies : (std::vector<double>) energies of the band matrix
		***************************************************/
		if(verbose == -1) std::cout<<"double bandLanczosAlgorithm(...) called"<<std::endl;
		//char sArrNature = (sArr->getLength('R') > 0 ? 'R':'K');
		//sArr->showAllStates();

		//Number of elements in array_bk
		int num_of_v = n_bk;
		product_c_omega->reserve(n_bk*5);
		std::vector<double> temp(n_bk* *nIter, 0);
		*product_c_omega = temp;
		//MinEnergy
		double previous_energy = 1000;
		std::vector<double> energies; 

		//Indexes of deflation
		std::vector<int> index_array; index_array.reserve(n_bk);
		//(1) Orthogonal basis
		double* zero = new double[len_bk]();
		std::vector<double> bk(vk->begin(), vk->end());
		
		//(2) number of maximum deflation
		int pc = n_bk;
		int M0 = 2 * n_bk +1;

		for (int i = pc; i < M0; i++)
			vk->insert(vk->end(), zero, zero + len_bk);

		int iterations = *nIter;

		//Elements for the new matrixes
		double* t_jpc = new double[iterations * iterations]();
		double* s_jpc = new double[iterations * iterations]();
		
		//Nature of treated states (used for application on the matrix)
		//bool same = sArr->allSameNature();
		//if(!same){sArrNature = 'M';}
		int j = 0;
		for (j = 0; j < iterations; j++){
			if(j%10==0 && verbose > 4) std::cout<< "Band Lanczos current iteration :"<< j << std::endl;	
			if(verbose > 99){
				for (int i = 0; i < M0; i++){
					double nn = cblas_dnrm2(len_bk, vk->data() + i * len_bk, 1);
					std::cout<<"vec["<<i<<"] = "<<toStringP_Q(nn)<<std::endl;
				}
			}
			//(3) Norm of the v_j vector
			double v_norm = cblas_dnrm2(len_bk, vk->data() + (j%M0) * len_bk, 1);
			//std::cout<<"NORM:"<<v_norm<<std::endl;
			//(4) Is the v_j vector negligeable
			if (v_norm <= dtol) {
				if (verbose > 9) std::cout<<"DELFLATION"<<std::endl;
				//Add index the deflated array
				if (j - pc >= 0) index_array.push_back(j - pc);//(a)
				//Erase current vector cause negligeable
				pc--;//(b)
				for (int q = 0; q < pc; q++) {
					std::copy(vk->begin() + ((j+1+q)%M0) * len_bk,vk->begin() + ((j+1+q)%M0 +1) * len_bk,vk->begin() + ((j+q)%M0) * len_bk);
				}
				std::copy(zero, zero + len_bk, vk->begin() + ((j+pc)%M0) * len_bk);
				num_of_v--;
				if (pc == 0) break;
				j--;
				continue;//(d)
			}
			//(5) Normalize v_j
			double t_m1 = 1 / v_norm;
			cblas_dscal(len_bk, t_m1, vk->data() + (j%M0) * len_bk, 1);
		
			//Add terms to t matrix 
			if (j >= pc) {t_jpc[j * iterations + j - pc] = v_norm;}

			//Qmatrix product requirements <phi|c_mu|Omega>
			double dot_product;
			for (uInt k = 0; k < n_bk; k++) {
				dot_product = cblas_ddot(len_bk, bk.data() + k * len_bk, 1, vk->data() + (j%M0) * len_bk, 1);
				
				(*product_c_omega)[k * *nIter + j] = dot_product;
			}

			//(6) Makes all the next vectors orthogonal to vj
			for (int k = j + 1; k < j + pc; k++) {
				//Dot product between vj and vk
				double vjvk;
				vjvk = cblas_ddot(len_bk, vk->data() + (j%M0) * len_bk, 1, vk->data() + (k%M0) * len_bk, 1);
				
				//Makes orthogonality
				double a = - vjvk;
				cblas_daxpy(len_bk, a, vk->data() + len_bk * (j%M0), 1, vk->data() + len_bk * (k%M0), 1);

				//Adding to the new element matrix
				if (k >= pc) {t_jpc[j * iterations + k - pc] = vjvk;}
			}

			//(7) Projection of the matrix
			//Applies matrix according to the States used
			std::copy(zero, zero + len_bk, vk->data() + (j + pc)%M0 * len_bk);
			sArr->H(vk->data() + ((j + pc)%M0) * len_bk, vk->data() + (j%M0) * len_bk);
			//std::cout<<"FROM"<<std::endl;
			//printVector(vk->data() + (j%M0)*len_bk, len_bk);
			//std::cout<<"APPLICATION OF H"<<std::endl;
			//printVector(vk->data() + ((j+pc)%M0)*len_bk,len_bk);
			//(8)Make sure that the new vector is othogonal to the previous ones
			int k0 = 0;
			if (k0 < j - pc) {k0 = j-pc;}
			for (int k = k0; k < j; k++){
				//Makes t_jpc hermitian
				t_jpc[k * iterations + j] = t_jpc[j * iterations + k];
				double a = -t_jpc[k * iterations + j];
				cblas_daxpy(len_bk, a, vk->data() + len_bk * (k%M0), 1, vk->data() + len_bk * ((j + pc)%M0), 1);
			}
		
			//(9) Removes from the new vector created in (9), the removed indexes and the current vector
			std::sort(index_array.begin(),index_array.end());

			////Deflated indexes
			for (unsigned long k = 0; k < index_array.size(); k++) {
				if(index_array.at(k) != j) continue;

				double dot_product;
				dot_product = cblas_ddot(len_bk, vk->data() + (index_array.at(k)%M0) * len_bk, 1, vk->data() + ((j + pc)%M0) * len_bk, 1);
				t_jpc[index_array.at(k) * iterations + j] = dot_product;

				double a = -t_jpc[index_array.at(k) * iterations + j];
				cblas_daxpy(len_bk, a, vk->data() + len_bk * (index_array.at(k)%M0), 1, vk->data() + len_bk * ((j + pc)%M0), 1);
			}

			////Diag element t(j,j)
			double VkVjpc, temp_minus;
			VkVjpc = cblas_ddot(len_bk, vk->data() + (j%M0) * len_bk, 1, vk->data() + ((j + pc)%M0) * len_bk, 1);
			t_jpc[j *iterations +j] = VkVjpc;
			
			temp_minus = -VkVjpc;
			cblas_daxpy(len_bk, temp_minus, vk->data() + len_bk * (j%M0), 1, vk->data() + len_bk * ((j + pc)%M0), 1);

			//(10) Manages Deflation
			for (unsigned long k = 0; k < index_array.size(); k++){
				s_jpc[j * iterations + index_array.at(k)] = t_jpc[index_array.at(k) * iterations + j]; 
			}
			if ((j+1)%n_bk == 0 && j>= ((int)n_bk-1)) {
				int jj = j+1;
				//(11) Creates the T_j matrix to solve
				double* T_jPr = new double[jj * jj]();
				for (int i = 0; i < jj; i++) {
					for (int l = i; l < jj; l++) {
						T_jPr[i * jj + l] = t_jpc[i * iterations + l] + s_jpc[i * iterations + l];
						if (i!=l) T_jPr[l * jj + i] = t_jpc[l * iterations + i] + s_jpc[l * iterations + i];
					}
				}
				//std::cout<<"T MATRIX"<<std::endl;
				//printMatrix(T_jPr,jj,jj,7,3);
				
				//(12) Solve T_j to check for convergence only eigen values
				double* eigen_values = new double[jj];
				//Tools for zheev
				char jobs = 'N', uplo='U';
				int lwork = (jj)*(jj+1);
				double* work = new double[lwork];
				double* rwork = new double[lwork];
				int info;
				dsyev_(&jobs,&uplo,&jj,T_jPr,&jj,eigen_values,work,&lwork,&info,1,1);
				delete[] T_jPr;
				//Delete zheev tools
				delete[] work; delete[] rwork;
				
				//std::cout<<"ZHEGV DONE"<<std::endl;
				//printVector(eigenValues,jj);
				double current_energy = eigen_values[0];	
				//Checks if the lowest eigen value has converged
				//std::cout<<"DELTA ENERGY:"<<currentEnergy - previousEnergy<<std::endl;
				//std::cout<<"CURR ENERGY:"<<toStringP_Q(currentEnergy)<<"\tPREV"<<toStringP_Q(previousEnergy)<<std::endl;
				if (abs(current_energy - previous_energy) < dtol || ((len_bk - jj) < n_bk)){
					delete[] eigen_values;
					break;
				}
				previous_energy = current_energy;
				delete[] eigen_values;  
			}//END OF IF
		}//End of For

		//Finds eigen vectors
		int jj = j;
		double* T_jPr = new double[jj * jj]();
		///Create the T_j matrix
		for (int i = 0; i < jj; i++) {
			for (int l = i; l < jj; l++) {
				T_jPr[i * jj + l] = t_jpc[i * iterations + l] + s_jpc[i * iterations + l];
				if (i!=l) T_jPr[l * jj + i] = t_jpc[l * iterations + i] + s_jpc[l * iterations + i];
			}
		}
		double* eigen_values = new double[jj];
		
		//printMatrix(T_jPr,jj,jj,5,3);
		///Tools for zheev
		char jobs = 'V', uplo='U';
		int lwork = (jj)*(jj+1);
		double* work = new double[lwork];
		double* rwork = new double[lwork];
		int info;

		dsyev_(&jobs,&uplo,&jj,T_jPr,&jj,eigen_values,work,&lwork,&info,1,1);
		///Delete tools for zheev
		delete[] work; delete[] rwork;
		//std::cout<<"ZHEGV DONE"<<std::endl;
		//printVector(eigenValues,jj);
		
		///Put the energies and the eigen vectors in vector
		energies = std::vector<double>(eigen_values,eigen_values + jj);
		*sub_space_vectors = std::vector<double>(T_jPr, T_jPr + jj * jj);

		if(verbose > 4) std::cout<<"Band Lanczos number of iteration until convergence : "<<j+1<<std::endl;
		for (int i = n_bk - 1; i >= 0; i--) {
			product_c_omega->erase(product_c_omega->begin() + jj + i * *nIter, product_c_omega->begin() + *nIter * (i + 1));
		}
		*nIter = jj;
		delete[] t_jpc; delete[] s_jpc;
		delete[] eigen_values; delete[] T_jPr;
		delete[] zero;
		return energies;
	}

	std::vector<double> bandLanczosAlgorithm_GE(std::vector<double>* vk, uInt n_bk, unsigned long len_bk, StatesArrType* sArr, uInt* nIter, std::vector<double>* subSpace_vectors, std::vector<double>* productCOmega, double dtol = 10e-10){
		std::vector<double> a;
		return a;
	}

	double fundEnergy(std::vector<double>* fundState, StatesArrType* states, int* deg){
		/***************************************************
		* Finds the fundamental energy by repeating the Lanczos algorithm until the minimum value has converged on a value
		*
		* Parameters
		* ----------
		* fundState : (std::vector<double>*) fundamental state of the given Hamiltonian matrix
		* states	: (StatesArrType*) array of states of the subSpace
		* deg		: (int*) degeneracy of the fundamental
		*
		* Returns
		* -------
		* fundEnergy: (double) fundamental energy 
		****************************************************/
		if(verbose == -1) std::cout<<"double fundEnergy(...) called"<<std::endl;
		double fundEnergy;
		int rows = states->getLength();

		if (rows > LANCZOS_SIZE) {
			fundEnergy = lanczosAlgorithm(fundState,states,deg);
		}
		else { //Will do the same as above put simpler because the matrix is much smaller
			double* H = new double[rows*rows]();
			//bool sameNat = states->allSameNature();

			states->matrixCreation(H);
			//printMatrix(H,rows,rows,3,3);

			char jobs = 'V', uplo='U';
			double* eigenValues = new double[rows];
			int lwork = rows*(rows+1);
			double* work = new double[lwork];
			double* rwork = new double[lwork];
			int info;
			dsyev_(&jobs,&uplo,&rows,H,&rows,eigenValues,work,&lwork,&info,1,1);

			fundEnergy = eigenValues[0];
			delete[] work; delete[] rwork;

			//Check degeneracy
			*deg = degFundamentalCheck(eigenValues,rows);
			delete[] eigenValues;

			if (*deg > 1) fundState->resize(rows*(*deg));
			//Stores the fundamental vector and if needed the degenerated ones too
			for (int j = 0; j < *deg; j++) {
				for(int i = 0; i<rows; i++){
					fundState->at(i+j*rows) = H[i+rows*j];
				}
			}
			//conjugateVector(fundState,rows);
			delete[] H;
		}

		return fundEnergy;
	}
};

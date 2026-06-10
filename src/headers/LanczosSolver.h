#pragma once

#include "ModulesStates/StatesR_T.h"
#include "ModulesStates/StatesK_T.h"
#include "ModulesStates/StatesR_H.h"
#include "basicFunctions.h"

int ONE = 1;
double ALPHA_D = 1;
double BETA_D = 0;
std::complex<double> ALPHA_C = 1;
std::complex<double> BETA_C = 0;

//GENERIC TEMPLATE
template<class T, class StatesArrType> class LanczosSolver;

//DOUBLE TEMPLATE
template<class StatesArrType> class LanczosSolver<double,StatesArrType>{
	public:
	void lanczos_energy(
        std::vector<double>* fundState_lanczos_basis, double* init_vector,
        StatesArrType* sArr, std::vector<double>* alpha,
        std::vector<double>* beta, double* fund_energy, int* iter, int* deg,
        double epsilon=10e-12) {
		/*******************************************************
		* Computes the alpha/beta of the tridiagonal matrix and find the
        * converged fund energy and fund vector in the reduced space.
		*
		* Parameters
		* ----------
		* fund_state_lanczos_basis	: (std::vector<double>*) fundamental eigen
        *                                           vector in the Lanczos basis
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
        //Number of states
		sType size = sArr->get_length();
		//Sets size of alpha and beta
		alpha->clear();	alpha->reserve(*iter);
		beta->clear();	beta->reserve(*iter);

		//Two main vectors
		std::vector<double> r(init_vector, init_vector + size);
		std::vector<double> q(size);

		//Energies to converge
		double prev_iter_energy = 10e10;
		double energy = 10e5;

		sType current_iteration = 1;
		bool converged = false;

		while (!converged) {
            //Normalization
			if (beta->size()) {
                double beta_m1 = 1/beta->back();
                double minus_beta = -beta->back();
			    dscal_(&size, &beta_m1, r.data(), &ONE);
			    dscal_(&size, &minus_beta, q.data(), &ONE);
			}

			std::vector<double> H_tmp(size);
			//Applies the vector r on the matrix H and stores it in H_tmp
			sArr->H(H_tmp.data(), r.data());

			daxpy_(&size, &ALPHA_D, H_tmp.data(), &ONE, q.data(), &ONE);//q = q + H*r

			//swap q <-> r
			dswap_(&size, q.data(), &ONE, r.data(), &ONE);

			double dot_product = ddot_(&size, q.data(), &ONE, r.data(), &ONE);
			alpha->push_back(dot_product);

            //r = r - q*alpha
            double minus_alpha = -alpha->back();
			daxpy_(&size, &minus_alpha, q.data(), &ONE, r.data(), &ONE);

			beta->push_back(dnrm2_(&size, r.data(), &ONE));

			//Arrays for tridiag solve
			double* arr_a = new double[alpha->size()];
			double* arr_b = new double[beta->size()];
			std::copy(alpha->begin(), alpha->end(), arr_a);
			std::copy(beta->begin(), beta->end(), arr_b);

			//Parameters for solver
			char jobs = 'V';
			int n = current_iteration;
			int info;
			double* vecs = new double[n*n];
			double* work = new double[2*n];

			//Solving Energy
			dstev_(&jobs, &n, arr_a, arr_b, vecs, &n, work, &info);
			//Fund energy

			prev_iter_energy = energy;
			energy = arr_a[0];


			if (abs(prev_iter_energy - energy) < epsilon
                && current_iteration > 3) {
				converged = true;
                //Check degeneracy
				*deg = deg_fundamental_check(arr_a, n);
				fundState_lanczos_basis->clear();
				*fundState_lanczos_basis = std::vector<double>(vecs,
                                                              vecs + n*(*deg));
			}
            else if (current_iteration == size) {
				converged = true;
                //Check degeneracy
				*deg = deg_fundamental_check(arr_a, n);
				fundState_lanczos_basis->clear();
				*fundState_lanczos_basis = std::vector<double>(vecs,
                                                              vecs + n*(*deg));
			}

			delete[] vecs; delete[] work; delete[] arr_a; delete[] arr_b;

			if (current_iteration%10 == 0 && verbose > 4) {
                printf("\rLanczos energy iteration : %4ld\tdE = %1.5e",current_iteration,abs(prev_iter_energy - energy));
                fflush(stdout);
			}
			current_iteration++;

		}
		if (verbose > 4) {
			std::cout << "\nLanczos iteration used : "
                << current_iteration << std::endl;
		}
		*iter = current_iteration;
		*fund_energy = energy;
	}
	void lanczos_vectors(
        std::vector<double>* fundState_lanczos_basis, double* gs,
        StatesArrType* sArr, std::vector<double>* alpha,
        std::vector<double>* beta, int* deg) {
		/*******************************************************
		* Convert the fund vector in the reduced space to the original space.
		*
		* Parameters
		* ----------
		* fund_state_lanczos_basis	: (std::vector<double>*) fundamental eigen
        *                                           vector in the Lanczos basis
		* gs						: (double) initial vector |phi_0> and on
        *                                                   out groundstate
		* sArr						: (StatesArrType) array of states object
		* alpha						: (std::vector<double>*) alphas
		* beta						: (std::vector<double>*) betas
		* deg						: (int*) degeneracy of the fund vector
		*
		* Returns
		* -------
		* NONE
		*******************************************************/
		sType size = sArr->get_length();
		int size_proj = alpha->size();

		std::vector<double> r(gs, gs + size);
		std::vector<double> q(size);

		for (int d = 0; d < *deg; d++) {
            double scal = fundState_lanczos_basis->at(size_proj*d);
            dscal_(&size, &scal, gs + d*size, &ONE);
        }


		for (unsigned int j = 1; j < fundState_lanczos_basis->size(); j++) {
			std::vector<double> H_tmp(size);
			sArr->H(H_tmp.data(), r.data());

			daxpy_(&size, &ALPHA_D, H_tmp.data(), &ONE, q.data(), &ONE);	//q = q + H*r

            //r = r - q*alpha
            double minus_alpha = -alpha->at(j-1);
			daxpy_(&size, &minus_alpha, r.data(), &ONE, q.data(), &ONE);
			for (unsigned int i = 0; i < r.size(); i++) {
				double tmp = r[i];
				r[i] = q[i]/beta->at(j-1);
				q[i] = -beta->at(j-1)*tmp;
			}

			for (int d = 0; d < *deg; d++) {
				double scal = fundState_lanczos_basis->at(j + size_proj*d);
                //r = r - q*alpha
				daxpy_(&size, &scal, r.data(), &ONE, gs + d*size, &ONE);
			}
			if (j%10 == 0 && verbose > 4) {

                print_iteration(j,"Lanczos vector iteration :");
			}
		}
	}
	double lanczos_algorithm(
        std::vector<double>* fundState, StatesArrType* sArr, int* deg,
        double epsilon = 10e-12) {
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
		if(verbose == -1) std::cout
            << "double lanczosAlgorithm(double...) called"<<std::endl;

		double fundEnergy;
		std::vector<double> alpha, beta;
		std::vector<double> fundState_lanczosBasis;
		int nIterations = 1000;

		// Random Initial Vector
		initial_vector(sArr->get_length(), fundState->data());

		lanczos_energy(&fundState_lanczosBasis, fundState->data(), sArr, &alpha,
                &beta, &fundEnergy, &nIterations, deg, epsilon);
		if (verbose > 9) std::cout<<"DEGENERACY:"<<*deg<<std::endl;

		//Increase the size of the fundState according to the degeneracy
		for (int i = 1; i < *deg; i++) {
			fundState->insert(fundState->end(), fundState->begin(),
                     fundState->begin() + sArr->get_length());
		}
		lanczos_vectors(&fundState_lanczosBasis, fundState->data(), sArr,
                 &alpha, &beta, deg);

		return fundEnergy;
	}


	std::vector<double> band_lanczos_algorithm(
        std::vector<double>* vk, uInt n_bk, sType len_bk,
        StatesArrType* sArr, uInt* nIter,
        std::vector<double>* sub_space_vectors,
        std::vector<double>* product_c_omega, double dtol = 10e-10){
		/*************************************************
		* Band lanczos algorithm
		*
		* Parameters
		* ----------
		* vk				: (std::vector<double>*) The n_bk initial vectors
		* n_bk				: (uInt) number of band vectors
		* len_bk			: (sType) length of the vectors
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
		if(verbose == -1) std::cout
            << "double bandLanczosAlgorithm(...) called"<<std::endl;

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

		int j = 0;
		for (j = 0; j < iterations; j++){
			if(j%10 == 0 && verbose > 4)
                print_iteration(j,"Band Lanczos iteration :");
			if(verbose > 99){
				for (int i = 0; i < M0; i++){
					double nn = dnrm2_(&len_bk, vk->data() + i*len_bk, &ONE);
					std::cout << "vec[" << i << "] = " << to_string_pq(nn)
                        << std::endl;
				}
			}
			//(3) Norm of the v_j vector
			double v_norm = dnrm2_(&len_bk, vk->data() + (j%M0)*len_bk, &ONE);
			//(4) Is the v_j vector negligeable
			if (v_norm <= dtol) {
				if (verbose > 9) std::cout << "DELFLATION" << std::endl;
				//Add index the deflated array
				if (j - pc >= 0) index_array.push_back(j - pc);//(a)
				//Erase current vector cause negligeable
				pc--;//(b)
				for (int q = 0; q < pc; q++) {
					std::copy(
                        vk->begin() + ((j+1+q)%M0) * len_bk,
                        vk->begin() + ((j+1+q)%M0 +1) * len_bk,
                        vk->begin() + ((j+q)%M0) * len_bk);
				}
				std::copy(zero, zero + len_bk,
                                        vk->begin() + ((j + pc)%M0) * len_bk);
				num_of_v--;
				if (pc == 0) break;
				j--;
				continue;//(d)
			}
			//(5) Normalize v_j
			double t_m1 = 1 / v_norm;
			dscal_(&len_bk, &t_m1, vk->data() + (j%M0) * len_bk, &ONE);

			//Add terms to t matrix
			if (j >= pc) {t_jpc[j * iterations + j - pc] = v_norm;}

			//Qmatrix product requirements <phi|c_mu|Omega>
			double dot_product;
			for (uInt k = 0; k < n_bk; k++) {
				dot_product = ddot_(&len_bk, bk.data() + k * len_bk, &ONE,
                                         vk->data() + (j%M0) * len_bk, &ONE);

				(*product_c_omega)[k * *nIter + j] = dot_product;
			}

			//(6) Makes all the next vectors orthogonal to vj
			for (int k = j + 1; k < j + pc; k++) {
				//Dot product between vj and vk
				double vjvk;
				vjvk = ddot_(&len_bk, vk->data() + (j%M0) * len_bk, &ONE,
                                  vk->data() + (k%M0) * len_bk, &ONE);

				//Makes orthogonality
				double a = - vjvk;
				daxpy_(&len_bk, &a, vk->data() + len_bk * (j%M0), &ONE,
                      vk->data() + len_bk * (k%M0), &ONE);

				//Adding to the new element matrix
				if (k >= pc) {t_jpc[j * iterations + k - pc] = vjvk;}
			}

			//(7) Projection of the matrix
			//Applies matrix according to the States used
			std::copy(zero, zero + len_bk, vk->data() + (j + pc)%M0 * len_bk);
			sArr->H(vk->data() + ((j + pc)%M0) * len_bk,
                    vk->data() + (j%M0) * len_bk);

			//(8)Make sure that the new vector is othogonal to the previous ones
			int k0 = 0;
			if (k0 < j - pc) {k0 = j-pc;}
			for (int k = k0; k < j; k++){
				//Makes t_jpc hermitian
				t_jpc[k * iterations + j] = t_jpc[j * iterations + k];
				double a = -t_jpc[k * iterations + j];
			    daxpy_(&len_bk, &a, vk->data() + len_bk * (k%M0), &ONE,
                      vk->data() + len_bk * ((j + pc)%M0), &ONE);
			}

			//(9) Removes from the new vector created in (9),
            //    the removed indexes and the current vector
			std::sort(index_array.begin(), index_array.end());

			////Deflated indexes
			for (unsigned long k = 0; k < index_array.size(); k++) {
				if(index_array.at(k) != j) continue;

				double dot_product;
				dot_product = ddot_(
                    &len_bk, vk->data() + (index_array.at(k)%M0) * len_bk, &ONE,
                    vk->data() + ((j + pc)%M0) * len_bk, &ONE);
				t_jpc[index_array.at(k) * iterations + j] = dot_product;

				double a = -t_jpc[index_array.at(k) * iterations + j];
				daxpy_(&len_bk, &a, vk->data() + len_bk * (index_array.at(k)%M0),
                      &ONE, vk->data() + len_bk * ((j + pc)%M0), &ONE);
			}

			////Diag element t(j,j)
			double VkVjpc, temp_minus;
			VkVjpc = ddot_(&len_bk, vk->data() + (j%M0) * len_bk, &ONE,
                                vk->data() + ((j + pc)%M0) * len_bk, &ONE);
			t_jpc[j *iterations +j] = VkVjpc;

			temp_minus = -VkVjpc;
			daxpy_(&len_bk, &temp_minus, vk->data() + len_bk * (j%M0), &ONE,
                  vk->data() + len_bk * ((j + pc)%M0), &ONE);

			//(10) Manages Deflation
			for (unsigned long k = 0; k < index_array.size(); k++) {
                double temp = t_jpc[index_array.at(k) * iterations + j];
				s_jpc[j * iterations + index_array.at(k)] = temp;
			}
			if ((j+1)%n_bk == 0 && j >= ((int)n_bk-1)) {
				int jj = j+1;
				//(11) Creates the T_j matrix to solve
				double* T_jPr = new double[jj * jj]();
				for (int i = 0; i < jj; i++) {
					for (int l = i; l < jj; l++) {
						T_jPr[i * jj + l] = t_jpc[i * iterations + l]
                                            + s_jpc[i * iterations + l];
                        //Symetric matrix
						if (i != l) {
                            T_jPr[l * jj + i] = t_jpc[l * iterations + i]
                                + s_jpc[l * iterations + i];
                        }
					}
				}

				//(12) Solve T_j to check for convergence only eigen values
				double* eigen_values = new double[jj];
				//Tools for dsyev
				char jobs = 'N', uplo='U';
				int lwork = (jj)*(jj+1);
				double* work = new double[lwork];
				double* rwork = new double[lwork];
				int info;

				dsyev_(&jobs, &uplo, &jj, T_jPr, &jj, eigen_values, work,
                       &lwork, &info);

				delete[] T_jPr;
				//Delete dsyev tools
				delete[] work; delete[] rwork;

				double current_energy = eigen_values[0];

				//Checks if the lowest eigen value has converged
				if (abs(current_energy - previous_energy) < dtol
                        || (((int)len_bk - jj) < (int)n_bk)) {
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
				T_jPr[i * jj + l] = t_jpc[i * iterations + l]
                                    + s_jpc[i * iterations + l];
				if (i != l) {
                    T_jPr[l * jj + i] = t_jpc[l * iterations + i]
                                        + s_jpc[l * iterations + i];
                }
			}
		}

		double* eigen_values = new double[jj];

		///Tools for dsyev
		char jobs = 'V', uplo='U';
		int lwork = (jj)*(jj+1);
		double* work = new double[lwork];
		double* rwork = new double[lwork];
		int info;

		dsyev_(&jobs,&uplo,&jj,T_jPr,&jj,eigen_values,work,&lwork,&info);
		///Delete tools for dsyev
		delete[] work; delete[] rwork;

		///Put the energies and the eigen vectors in vector
		energies = std::vector<double>(eigen_values, eigen_values + jj);
		*sub_space_vectors = std::vector<double>(T_jPr, T_jPr + jj * jj);

		if(verbose > 4) std::cout << "\nBand Lanczos number of iteration until"
                                  << " convergence : "<< j+1 << std::endl;
		for (int i = n_bk - 1; i >= 0; i--) {
			product_c_omega->erase(product_c_omega->begin() + jj + i * *nIter,
                          product_c_omega->begin() + *nIter * (i + 1));
		}
		*nIter = jj;
		delete[] t_jpc; delete[] s_jpc;
		delete[] eigen_values; delete[] T_jPr;
		delete[] zero;
		return energies;
	}

	double fund_energy(std::vector<double>* fund_state, StatesArrType* states, int* deg){
		/***************************************************
		* Finds the fundamental energy by repeating the Lanczos algorithm until the minimum value has converged on a value
		*
		* Parameters
		* ----------
		* fund_state : (std::vector<double>*) fundamental state of the given Hamiltonian matrix
		* states	: (StatesArrType*) array of states of the subSpace
		* deg		: (int*) degeneracy of the fundamental
		*
		* Returns
		* -------
		* fund_energy: (double) fundamental energy
		****************************************************/
		if(verbose == -1) std::cout << "double fundEnergy(...) called\n";
		double fund_energy;
		int rows = states->get_length();

		if (rows > LANCZOS_SIZE) {
			fund_energy = lanczos_algorithm(fund_state, states, deg);
		}
		else {
            //Will do the same as above put with a smaller matrix
			double* H = new double[rows*rows]();

			states->matrix_creation(H);

			char jobs = 'V', uplo='U';
			double* eigen_values = new double[rows];
			int lwork = rows*(rows+1);
			double* work = new double[lwork];
			double* rwork = new double[lwork];
			int info;
			dsyev_(&jobs, &uplo, &rows, H, &rows, eigen_values, work, &lwork,
                   &info);

			fund_energy = eigen_values[0];
			delete[] work; delete[] rwork;

			//Check degeneracy
			*deg = deg_fundamental_check(eigen_values, rows);
			delete[] eigen_values;

			if (*deg > 1) fund_state->resize(rows*(*deg));
			//Stores the fundamental vector and
            //if needed the degenerated ones too
			for (int j = 0; j < *deg; j++) {
				for(int i = 0; i<rows; i++){
					fund_state->at(i+j*rows) = H[i+rows*j];
				}
			}
			delete[] H;
		}

		return fund_energy;
	}
};
//
//COMPLEX DOUBLE TEMPLATE
template<class StatesArrType> class LanczosSolver<std::complex<double>,StatesArrType>{
	public:
	void lanczos_energy(std::vector<double>* fundState_lanczosBasis, std::complex<double>* initVector, StatesArrType* sArr, std::vector<double>* alpha, std::vector<double>* beta, double* fundEnergy, int* iter, int* deg, double epsilon=10e-12) {
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
        //Number of states
		sType size = sArr->get_length();
		//Sets size of alpha and beta
		alpha->clear();	alpha->reserve(*iter);
		beta->clear();	beta->reserve(*iter);

		//Two main vectors
		std::vector<std::complex<double>> r(initVector, initVector + size);
		std::vector<std::complex<double>> q(size);

		//Energies to converge
		double prevIterEnergy = 1000;
		double energy = 100;

		int currentIteration = 1;

		bool converged = false;

		while (!converged) {
			if(beta->size()) {
                double beta_m1 = 1/beta->back();
                double minus_beta = -beta->back();
				zdscal_(&size, &beta_m1, r.data(), &ONE);
				zdscal_(&size, &minus_beta, q.data(), &ONE);
			}

			std::vector<std::complex<double>> H_tmp(size);
		//Applies the vector r on the matrix H and stores it in H_tmp
			sArr->H(H_tmp.data(),r.data());

			std::complex<double> one = 1;
			zaxpy_(&size, &one, H_tmp.data(),&ONE,q.data(),&ONE);	//q = q + H*r
			//swap q <-> r
		    zswap_(&size, q.data(),&ONE,r.data(),&ONE);

			std::complex<double> dotProd = zdotc_(&size, q.data(), &ONE, r.data(), &ONE);
			alpha->push_back(dotProd.real());

			std::complex<double> alpha_back = -alpha->back();
			zaxpy_(&size, &alpha_back, q.data(), &ONE, r.data(), &ONE);	//r = r - q*alpha

			beta->push_back(dznrm2_(&size, r.data(), &ONE));

			//Arrays for tridiag solve
			double* arr_a = new double[alpha->size()];
			double* arr_b = new double[beta->size()];
			std::copy(alpha->begin(),alpha->end(),arr_a);
			std::copy(beta->begin(),beta->end(),arr_b);

			//Parameters for solver
			char jobs = 'V';
			int n = currentIteration;
			int info;
			double* vecs = new double[n*n];
			double* work = new double[2*n];

			//Solving Energy
			dstev_(&jobs,&n,arr_a,arr_b,vecs,&n,work,&info);
			//Fund energy

			prevIterEnergy = energy;
			energy = arr_a[0];
			double RITZ = std::abs(vecs[n-1]*beta->back());
            //std::cout<<"RITZ VALUE : "<< RITZ<<std::endl;
			if (abs(prevIterEnergy - energy) < epsilon && currentIteration > 3) {
				converged = true;
				*deg = deg_fundamental_check(arr_a,n);
				fundState_lanczosBasis->clear();
				*fundState_lanczosBasis = std::vector<double>(vecs, vecs + n*(*deg));
			}

			delete[] vecs; delete[] work; delete[]arr_a; delete[] arr_b;

			if(currentIteration%10==0 && verbose > 5){
				std::cout << "Lanczos current iteration :" << currentIteration << ":"<<abs(prevIterEnergy - energy)<< std::endl;

			}
			if(currentIteration == size) {
				break;
			}
			currentIteration++;

		}
		if(verbose > 5){
			std::cout << "Lanczos iteration used:" << currentIteration << std::endl;
		}
		*iter = currentIteration;
		*fundEnergy = energy;
	}

	void lanczos_vectors(std::vector<double>* fundState_lanczosBasis, std::complex<double>* gs, StatesArrType* sArr, std::vector<double>* alpha, std::vector<double>* beta, int* deg) {
		sType size = sArr->get_length();
		int size_proj = alpha->size();

		std::vector<std::complex<double>> r(gs, gs + size);
		std::vector<std::complex<double>> q(size);

		for (int d = 0; d < *deg; d++)
	        zdscal_(&size, &fundState_lanczosBasis->at(0+d*size_proj), gs+d*size, &ONE);

		for (unsigned int j = 1; j <fundState_lanczosBasis->size(); j++) {
			std::vector<std::complex<double>> H_tmp(size);
			sArr->H(H_tmp.data(),r.data());

			std::complex<double> scal = 1;
			zaxpy_(&size, &scal, H_tmp.data(), &ONE, q.data(), &ONE);	//q = q + H*r
			scal = -alpha->at(j-1);
			zaxpy_(&size, &scal, r.data(), &ONE, q.data(), &ONE);	//r = r - q*alpha
			for (unsigned int i = 0; i < r.size(); i++) {
				std::complex<double> tmp = r[i];
				r[i] = q[i]/beta->at(j-1);
				q[i] = -beta->at(j-1)*tmp;
			}

			for (int d = 0; d < *deg; d++) {
				std::complex<double> temp_fundState_lanczosBasis = fundState_lanczosBasis->at(j+size_proj*d);
				zaxpy_(&size, &temp_fundState_lanczosBasis, r.data(), &ONE, gs+d*size, &ONE);	//r = r - q*alpha
			}
		}
	}

	double lanczos_algorithm(std::vector<std::complex<double>>* fundState, StatesArrType* sArr, int* deg, double epsilon = 10e-12) {
		/*************************************************
		Redefines a given matrix with the Lanczos Algorithm without needing the Hamiltonian matrix

		Parameters
		----------
		fundState : (double) fundamental state of the system
		sArr : (StatesArr) States used in the subspace
		deg			: (int*) degeneracy counter
		epsilon : (double) Convergence acceptability

		Returns
		-------
		currentEnergy : (double) fundamental energy of the system
		***************************************************/
		if(verbose == -1) std::cout<<"double lanczosAlgorithm(complex...) called"<<std::endl;

		double fund_energy;
		std::vector<double> alpha, beta;
		std::vector<double> fundState_lanczosBasis;
		int nIterations = 1000;

		// Random Initial Vector
		initial_vector(sArr->get_length(),fundState->data());

		lanczos_energy(&fundState_lanczosBasis,fundState->data(),sArr,&alpha,&beta,&fund_energy,&nIterations,deg,epsilon);
		//Increase the size of the fundState according to the degeneracy
		for (int i = 1; i < *deg; i++) {
			fundState->insert(fundState->end(),fundState->begin(),fundState->begin()+sArr->get_length());
		}

		lanczos_vectors(&fundState_lanczosBasis,fundState->data(),sArr,&alpha,&beta,deg);

		return fund_energy;
	}

	std::vector<double> band_lanczos_algorithm(std::vector<std::complex<double>>* vk, uInt n_bk, unsigned long len_bk, StatesArrType* sArr, uInt* nIter, std::vector<std::complex<double>>* subSpace_vectors, std::vector<std::complex<double>>* productCOmega, double dtol = 10e-10){
		if(verbose == -1) std::cout<<"complex bandLanczosAlgorithm(...) called"<<std::endl;

		//char sArrNature = (sArr->get_length('R') > 0 ? 'R':'K');

		//Number of elements in array_bk
		int numOfV = n_bk;
		productCOmega->reserve(n_bk*5);
		std::vector<std::complex<double>> temp(n_bk* *nIter, 0);
		*productCOmega = temp;
		//MinEnergy
		double previousEnergy = 1000;
		std::vector<double> energies;

		//Indexes of deflation
		std::vector<int> index_array; index_array.reserve(n_bk);
		//(1) Orthogonal basis
		std::complex<double>* zero = new std::complex<double>[len_bk]();
		std::vector<std::complex<double>> bk(vk->begin(), vk->end());

		//(2) number of maximum deflation
		int pc = n_bk;
		int M0 = 2 * n_bk +1;

		for (int i = pc; i < M0; i++)
			vk->insert(vk->end(), zero, zero + len_bk);

		int iterations = *nIter;

		//Elements for the new matrixes
		std::complex<double>* t_jpc = new std::complex<double>[iterations * iterations]();
		std::complex<double>* s_jpc = new std::complex<double>[iterations * iterations]();

		//Nature of treated states (used for application on the matrix)
		//bool same = sArr->allSameNature();
		//if(!same){sArrNature = 'M';}
		int j = 0;
		for (j = 0; j < iterations; j++){
			if(j%10==0) if(verbose > 5) std::cout<< "Band Lanczos current iteration :"<< j << std::endl;
			if(verbose > 99){
				for (int i = 0; i < M0; i++){
					double nn = dznrm2_(&len_bk, vk->data() + i * len_bk, &ONE);
					std::cout<<"vec["<<i<<"] = "<<nn<<std::endl;
				}
			}
			//(3) Norm of the v_j vector
			double v_norm = dznrm2_(&len_bk, vk->data() + (j%M0) * len_bk, &ONE);
			//(4) Is the v_j vector negligeable
			if (v_norm <= dtol) {
				//Add index the deflated array
				if (j - pc >= 0) index_array.push_back(j - pc);//(a)
				//Erase current vector cause negligeable
				pc--;//(b)
				for (int q = 0; q < pc; q++) {
					std::copy(vk->begin() + ((j+1+q)%M0) * len_bk,vk->begin() + ((j+1+q)%M0+1) * len_bk,vk->begin() + ((j+q)%M0) * len_bk);
				}
					std::copy(zero, zero + len_bk, vk->begin() + ((j+pc)%M0) * len_bk);
				numOfV--;
				if (pc == 0) break;
				j--;
				continue;//(d)
			}
			//(5) Normalize v_j
			double t_m1 = 1 / v_norm;
			zdscal_(&len_bk, &t_m1, vk->data() + (j%M0) * len_bk, &ONE);

			//Add terms to t matrix
			if (j >= pc) {t_jpc[j * iterations + j - pc] = v_norm;}

			//Qmatrix product requirements <phi|c_mu|Omega>
			std::complex<double> dotProd;
			for (uInt k = 0; k < n_bk; k++) {
				dotProd = zdotc_(&len_bk, bk.data() + k * len_bk, &ONE, vk->data() + (j%M0) * len_bk, &ONE);
				(*productCOmega)[k * *nIter + j] = dotProd;
			}

			//(6) Makes all the next vectors orthogonal to vj
			for (int k = j + 1; k < j + pc; k++) {
				//Dot product between vj and vk
				std::complex<double> vjvk;
				vjvk = zdotc_(&len_bk, vk->data() + (j%M0) * len_bk, &ONE, vk->data() + (k%M0) * len_bk, &ONE);

				//Makes orthogonality
				std::complex<double> a = - vjvk;
				zaxpy_(&len_bk, &a, vk->data() + len_bk * (j%M0), &ONE, vk->data() + len_bk * (k%M0), &ONE);

				//Adding to the new element matrix
				if (k >= pc) {t_jpc[j * iterations + k - pc] = vjvk;}
			}

			//(7) Projection of the matrix
			//Applies matrix according to the States used
			std::copy(zero, zero + len_bk, vk->data() + (j + pc)%M0 * len_bk);
			sArr->H(vk->data() + ((j + pc)%M0) * len_bk, vk->data() + (j%M0) * len_bk);
			//std::cout<<"APPLICATION OF H"<<std::endl;
			//printVector(vk->data() + ((j+pc)%M0)*len_bk,len_bk);

			//(8)Make sure that the new vector is othogonal to the previous ones
			int k0 = 0;
			if (k0 < j - pc) {k0 = j-pc;}
			for (int k = k0; k < j; k++){
				//Makes t_jpc hermitian
				t_jpc[k * iterations + j] = conjugate(t_jpc[j * iterations + k]);
				std::complex<double> a = -t_jpc[k * iterations + j];
				zaxpy_(&len_bk, &a, vk->data() + len_bk * (k%M0), &ONE, vk->data() + len_bk * ((j + pc)%M0), &ONE);
			}

			//(9) Removes from the new vector created in (9), the removed indexes and the current vector
			std::sort(index_array.begin(),index_array.end());

			////Deflated indexes
			for (unsigned long k = 0; k < index_array.size(); k++) {
				if(index_array.at(k) != j) continue;

				std::complex<double> dot_product;
				dot_product = zdotc_(&len_bk, vk->data() + (index_array.at(k)%M0) * len_bk, &ONE, vk->data() + ((j + pc)%M0) * len_bk, &ONE);
				t_jpc[index_array.at(k) * iterations + j] = dot_product;

				std::complex<double> a = -t_jpc[index_array.at(k) * iterations + j];
			    zaxpy_(&len_bk, &a, vk->data() + len_bk * (index_array.at(k)%M0), &ONE, vk->data() + len_bk * ((j + pc)%M0), &ONE);
			}

			////Diag element t(j,j)
			std::complex<double> VkVjpc, tempMinus;
			VkVjpc = zdotc_(&len_bk, vk->data() + (j%M0) * len_bk, &ONE, vk->data() + ((j + pc)%M0) * len_bk, &ONE);
			t_jpc[j *iterations +j] = VkVjpc;

			tempMinus = -VkVjpc;
			zaxpy_(&len_bk, &tempMinus, vk->data() + len_bk * (j%M0), &ONE, vk->data() + len_bk * ((j + pc)%M0), &ONE);

			//(10) Manages Deflation
			for (unsigned long k = 0; k < index_array.size(); k++){
				s_jpc[j * iterations + index_array.at(k)] = conjugate(t_jpc[index_array.at(k) * iterations + j]);
			}
			if ((j+1)%n_bk == 0 && j>= ((int)n_bk-1)) {
				int jj = j+1;
				//(11) Creates the T_j matrix to solve
				std::complex<double>* T_jPr = new std::complex<double>[jj * jj]();
				for (int i = 0; i < jj; i++) {
					for (int l = i; l < jj; l++) {
						T_jPr[i * jj + l] = t_jpc[i * iterations + l] + s_jpc[i * iterations + l];
						if (i!=l) T_jPr[l * jj + i] = t_jpc[l * iterations + i] + s_jpc[l * iterations + i];
					}
				}

				//(12) Solve T_j to check for convergence only eigen values
				double* eigenValues = new double[jj];
				//Tools for zheev
				char jobs = 'N', uplo='U';
				int lwork = (jj)*(jj+1);
				std::complex<double>* work = new std::complex<double>[lwork];
				double* rwork = new double[lwork];
				int info;

				zheev_(&jobs, &uplo, &jj, T_jPr, &jj,eigenValues, work, &lwork, rwork, &info);
				delete[] T_jPr;
				//Delete zheev tools
				delete[] work; delete[] rwork;

				double currentEnergy = eigenValues[0];
				//Checks if the lowest eigen value has converged
				if (abs(currentEnergy - previousEnergy) < dtol || ((len_bk - jj) < n_bk)){
					delete[] eigenValues;
					break;
				}
				previousEnergy = currentEnergy;
				delete[] eigenValues;
			}//END OF IF
		}//End of For

		//Finds eigen vectors
		int jj = j;
		std::complex<double>* T_jPr = new std::complex<double>[jj * jj]();
		///Create the T_j matrix
		for (int i = 0; i < jj; i++) {
			for (int l = i; l < jj; l++) {
				T_jPr[i * jj + l] = t_jpc[i * iterations + l] + s_jpc[i * iterations + l];
				if (i!=l) T_jPr[l * jj + i] = t_jpc[l * iterations + i] + s_jpc[l * iterations + i];
			}
		}
		double* eigenValues = new double[jj];

		///Tools for zheev
		char jobs = 'V', uplo='U';
		int lwork = (jj)*(jj+1);
		std::complex<double>* work = new std::complex<double>[lwork];
		double* rwork = new double[lwork];
		int info;

		zheev_(&jobs, &uplo, &jj, T_jPr, &jj, eigenValues, work, &lwork, rwork, &info);
		///Delete tools for zheev
		delete[] work; delete[] rwork;

		///Put the energies and the eigen vectors in vector
		energies = std::vector<double>(eigenValues,eigenValues + jj);
		*subSpace_vectors = std::vector<std::complex<double>>(T_jPr, T_jPr + jj * jj);

		if(verbose > 4 && (j+1)%10==0) std::cout<<"Band Lanczos number of iteration until convergence : "<<j+1<<std::endl;
		for (int i = n_bk - 1; i >= 0; i--) {
			productCOmega->erase(productCOmega->begin() + jj + i * *nIter, productCOmega->begin() + *nIter * (i + 1));
		}
		*nIter = jj;
		delete[] t_jpc; delete[] s_jpc;
		delete[] eigenValues; delete[] T_jPr;
		delete[] zero;
		return energies;
	}


	double fund_energy(std::vector<std::complex<double>>* fundState, StatesArrType* states, int* deg) {
		/*****************************************************
		Finds the fundamental energy by repeating the Lanczos algorithm until the
		minimum value has converged on a value

		Parameters
		----------
		H: (std::complex<double> *) Hamiltonian matrix
		fundState : (std::complex<double> *) fundamental state of the given Hamiltonian matrix
		states : (StatesArr) states of the subSpace

		Returns
		-------
		fundEnergy : (double) fundamental energy of the given Hamiltonian matrix
		*******************************************************/
		if(verbose == -1) std::cout<<"complex fundEnergy(...) called"<<std::endl;
		double fundEnergy;
		int rows = states->get_length();

		if (rows > LANCZOS_SIZE) {
			//OLD
			//fundEnergy = lanczosAlgorithm(H, fundState, rows);
			//fundEnergy = lanczosAlgorithm_new2(fundState,states);
            fundEnergy = lanczos_algorithm(fundState,states,deg);
		}
		else { //Will do the same as above put simpler because the matrix is much smaller
			std::complex<double>* H = new std::complex<double>[rows*rows]();
			//bool sameNat = states->allSameNature();
			states->matrix_creation(H);

            print_matrix(H,rows,rows,2,2);


			char jobs = 'V', uplo='U';
			double* eigenValues = new double[rows];
			int lwork = rows*(rows+1);
			std::complex<double>* work = new std::complex<double>[lwork];
			double* rwork = new double[lwork];
			int info;
            zheev_(&jobs, &uplo, &rows, H, &rows, eigenValues,work, &lwork, rwork, &info);

			fundEnergy = eigenValues[0];
			delete[] work; delete[] rwork;

			*deg = deg_fundamental_check(eigenValues,rows);
			delete[] eigenValues;

			if (*deg > 1) fundState->resize(rows*(*deg));
			//Stores the fundamental vector and if needed the degenerated ones too
			for (int j = 0; j < *deg; j++) {
				for(int i = 0; i<rows; i++){
					fundState->at(i+j*rows) = H[i+rows*j];
				}
			}
			conjugate_vector(fundState->data(),rows*(*deg));
			delete[] H;
		}

		return fundEnergy;
	}

};


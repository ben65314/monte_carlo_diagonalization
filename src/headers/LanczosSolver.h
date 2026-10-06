#pragma once

#include "ModulesStates/StatesR_T.h"
#include "ModulesStates/StatesK_T.h"
#include "ModulesStates/StatesR_H.h"

//GENERIC TEMPLATE
//template<class T, class StatesArrType> class LanczosSolver;

//DOUBLE TEMPLATE
template<class D, class StatesArrType> class LanczosSolver{
	public:
	void lanczos_energy(
        std::vector<double>* fundState_lanczos_basis, D* init_vector,
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
		std::vector<D> r(init_vector, init_vector + size);
		std::vector<D> q(size);

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
			    scal_blas(&size, &beta_m1, r.data());
			    scal_blas(&size, &minus_beta, q.data());
			}

			std::vector<D> H_tmp(size);
			//Applies the vector r on the matrix H and stores it in H_tmp
			sArr->H(H_tmp.data(), r.data());

			axpy_blas(&size, H_tmp.data(), q.data());//q = q + H*r

			//swap q <-> r
			swap_blas(&size, q.data(), r.data());

            std::complex<double> dot_product = dot_blas(&size, q.data(), r.data());
			alpha->push_back(dot_product.real());

            //r = r - q*alpha
            D minus_alpha = -alpha->back();
			axpy_blas(&size, q.data(), r.data(), minus_alpha);

			beta->push_back(nrm2_blas(&size, r.data()));

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


			if (abs(prev_iter_energy - energy) < epsilon && current_iteration > 3) {
				converged = true;
                //Check degeneracy
				*deg = deg_fundamental_check(arr_a, n);
				fundState_lanczos_basis->clear();
				*fundState_lanczos_basis = std::vector<double>(vecs,
                                                              vecs + n*(*deg));
			} else if (current_iteration == size) {
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
            //if (current_iteration > 5) break;
		}

		if (verbose > 4) {
			std::cout << "\nLanczos iteration used : "
                << current_iteration << std::endl;
		}
		*iter = current_iteration;
		*fund_energy = energy;
	}
	void lanczos_vectors(
        std::vector<double>* fundState_lanczos_basis, D* gs,
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
		sType size_proj = alpha->size();

		std::vector<D> r(gs, gs + size);
		std::vector<D> q(size);

		for (int d = 0; d < *deg; d++) {
            double scal = fundState_lanczos_basis->at(size_proj*d);
            scal_blas(&size, &scal, gs + d*size);
        }


		for (sType j = 1; j < size_proj; j++) {
			std::vector<D> H_tmp(size);
			sArr->H(H_tmp.data(), r.data());

			axpy_blas(&size, H_tmp.data(), q.data());	//q = q + H*r

            //r = r - q*alpha
            D minus_alpha = -alpha->at(j-1);
			axpy_blas(&size, r.data(), q.data(), minus_alpha);
			for (unsigned int i = 0; i < r.size(); i++) {
				D tmp = r[i];
				r[i] = q[i]/beta->at(j-1);
				q[i] = -beta->at(j-1)*tmp;
			}

			for (int d = 0; d < *deg; d++) {
				D scal = fundState_lanczos_basis->at(j + size_proj*d);
                //r = r - q*alpha
				axpy_blas(&size, r.data(), gs + d*size, scal);
			}
			if (j%10 == 0 && verbose > 4) {

                print_iteration(j,"Lanczos vector iteration :");
			}
		}
	}
	double lanczos_algorithm(
        std::vector<D>* fundState, StatesArrType* sArr, int* deg,
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
        std::vector<D>* vk, uInt n_bk, sType len_bk,
        StatesArrType* sArr, uInt* nIter,
        std::vector<D>* sub_space_vectors,
        std::vector<D>* product_c_omega, double dtol = 10e-10){
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
		product_c_omega->reserve(n_bk*5);
		std::vector<D> temp(n_bk* *nIter, 0);
		*product_c_omega = temp;

		//MinEnergy
		double previous_energy = 1000;
		std::vector<double> energies;

		//Indexes of deflation
		std::vector<int> index_array; index_array.reserve(n_bk);
		//(1) Orthogonal basis
		double* zero = new double[len_bk]();
		std::vector<D> bk(vk->begin(), vk->end());

		//(2) number of maximum deflation
		int pc = n_bk;
		int M0 = 2 * n_bk +1;

		for (int i = pc; i < M0; i++)
			vk->insert(vk->end(), zero, zero + len_bk);

		int iterations = *nIter;

		//Elements for the new matrixes
		D* t_jpc = new D[iterations * iterations]();
		D* s_jpc = new D[iterations * iterations]();

		int j = 0;
		for (j = 0; j < iterations; j++){

			if(j%10 == 0 && verbose > 4)
                print_iteration(j,"Band Lanczos iteration :");
			if(verbose > 99){
				for (int i = 0; i < M0; i++){
					double nn = nrm2_blas(&len_bk, vk->data() + i*len_bk);
					std::cout << "vec[" << i << "] = " << to_string_pq(nn)
                        << std::endl;
				}
			}
			//(3) Norm of the v_j vector
			double v_norm = nrm2_blas(&len_bk, vk->data() + (j%M0)*len_bk);
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
				if (pc == 0) break;
				j--;
				continue;//(d)
			}
			//(5) Normalize v_j
			double t_m1 = 1 / v_norm;
			scal_blas(&len_bk, &t_m1, vk->data() + (j%M0) * len_bk);

			//Add terms to t matrix
			if (j >= pc) {t_jpc[j * iterations + j - pc] = v_norm;}

			//Qmatrix product requirements <phi|c_mu|Omega>
			D dot_product;
			for (uInt k = 0; k < n_bk; k++) {
				dot_product = dot_blas(&len_bk, bk.data() + k * len_bk,
                                         vk->data() + (j%M0) * len_bk);

				(*product_c_omega)[k * *nIter + j] = dot_product;
			}

			//(6) Makes all the next vectors orthogonal to vj
			for (int k = j + 1; k < j + pc; k++) {
				//Dot product between vj and vk
				D vjvk;
				vjvk = dot_blas(&len_bk, vk->data() + (j%M0) * len_bk,
                                  vk->data() + (k%M0) * len_bk);

				//Makes orthogonality
				D a = - vjvk;
				axpy_blas(&len_bk, vk->data() + len_bk * (j%M0),
                      vk->data() + len_bk * (k%M0), a);

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
				D a = -t_jpc[k * iterations + j];
			    axpy_blas(&len_bk, vk->data() + len_bk * (k%M0),
                      vk->data() + len_bk * ((j + pc)%M0), a);
			}

			//(9) Removes from the new vector created in (9),
            //    the removed indexes and the current vector
			std::sort(index_array.begin(), index_array.end());

			////Deflated indexes
			for (unsigned long k = 0; k < index_array.size(); k++) {
				if(index_array.at(k) != j) continue;

				D dot_product;
				dot_product = dot_blas(
                    &len_bk, vk->data() + (index_array.at(k)%M0) * len_bk,
                    vk->data() + ((j + pc)%M0) * len_bk);
				t_jpc[index_array.at(k) * iterations + j] = dot_product;

				D a = -t_jpc[index_array.at(k) * iterations + j];
				axpy_blas(&len_bk, vk->data() + len_bk * (index_array.at(k)%M0),
                            vk->data() + len_bk * ((j + pc)%M0), a);
			}

			////Diag element t(j,j)
			D VkVjpc, temp_minus;
			VkVjpc = dot_blas(&len_bk, vk->data() + (j%M0) * len_bk,
                                vk->data() + ((j + pc)%M0) * len_bk);
			t_jpc[j *iterations +j] = VkVjpc;

			temp_minus = -VkVjpc;
			axpy_blas(&len_bk, vk->data() + len_bk * (j%M0),
                  vk->data() + len_bk * ((j + pc)%M0), temp_minus);

			//(10) Manages Deflation
			for (unsigned long k = 0; k < index_array.size(); k++) {
				s_jpc[j * iterations + index_array.at(k)] = conjugate(t_jpc[index_array.at(k) * iterations + j]);
			}
			if ((j+1)%n_bk == 0 && j >= ((int)n_bk-1)) {
				int jj = j+1;
				//(11) Creates the T_j matrix to solve
				D* T_jPr = new D[jj * jj]();
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
				D* work = new D[lwork];
				double* rwork = new double[lwork];

                sType row = jj;
                heev_lapack(jobs, uplo, &row, T_jPr, eigen_values);

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
		D* T_jPr = new D[jj * jj]();
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

        sType row = jj;
        heev_lapack(jobs, uplo, &row, T_jPr, eigen_values);
		///Delete tools for dsyev
		delete[] work; delete[] rwork;

		///Put the energies and the eigen vectors in vector
		energies = std::vector<double>(eigen_values, eigen_values + jj);
		*sub_space_vectors = std::vector<D>(T_jPr, T_jPr + jj * jj);

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

	double fund_energy(std::vector<D>* fund_state, StatesArrType* states, int* deg){
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
		sType rows = states->get_length();

		if (rows > LANCZOS_SIZE) {
			fund_energy = lanczos_algorithm(fund_state, states, deg);
		}
		else {
            //Will do the same as above put with a smaller matrix
			D* H = new D[rows*rows]();

			states->matrix_creation(H);

			char jobs = 'V', uplo='U';
			double* eigen_values = new double[rows];
			int lwork = rows*(rows+1);
			double* work = new double[lwork];
			double* rwork = new double[lwork];
            heev_lapack(jobs, uplo, &rows, H, eigen_values);

			fund_energy = eigen_values[0];
			delete[] work; delete[] rwork;

			//Check degeneracy
			*deg = deg_fundamental_check(eigen_values, rows);
			delete[] eigen_values;

			if (*deg > 1) fund_state->resize(rows*(*deg));
			//Stores the fundamental vector and
            //if needed the degenerated ones too
			for (int j = 0; j < *deg; j++) {
				for(sType i = 0; i < rows; i++){
					fund_state->at(i+j*rows) = H[i+rows*j];
				}
			}
			delete[] H;
		}

		return fund_energy;
	}
};


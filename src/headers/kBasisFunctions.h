#pragma once
#include "electronManipulationFunctions.h"
// k-basis
template <class T>void HuN(T state, std::vector<T>* proj_states, int sites) {
	/************************************************************
	* Generates the possible transition of the interaction hubbard hamiltonian in k base
	*
	* Parameters
	* ----------
	* state			: (T) initial state in the Fock formalisme
	* proj_states	: (std::vector<T>*) Array of the possible accessible states after the Hamiltonian
	* sites			: (int) Number of sites of the system
	*
	* Templates
	* ---------
	* T			: int, short, long, unsigned
	*
	* Returns
	* -------
	* NONE
	**************************************************************/
	//States after applied Hamiltonian
	T iter_up = 1UL << sites;
	//Checks for all k
	for (char k = sites - 1; k >= 0; k--) {
		//Check for all destruction sites
		T iter_down = 1UL;
		for (char l = sites - 1; l >= 0; l--) {
			T evolving_state = state;
			//T change = 0;

			if (((evolving_state & iter_up) != 0) && ((evolving_state & iter_down) != 0)) {

				evolving_state ^= (iter_up | iter_down);
				//change -= (iterUp + iterDown);
			}
			else {iter_down <<= 1;	continue;}
			for (unsigned char q = 0; q < sites; q++) {
				T final_state;
				int qk = q+k;
				T where_up;
				if (qk >= sites) {
					where_up = iter_up << (sites - q);
				}
				else where_up = iter_up >> (q);

				int ql = l - q;
				T where_down;
				if (ql < 0) {
					where_down = iter_down >> (sites - q);

				}
				else where_down = iter_down << q;


				if (((evolving_state | where_up) != evolving_state) && ((evolving_state | where_down) != evolving_state)) {
					final_state = evolving_state | (where_up | where_down);
					//change += (whereUp + whereDown);
				}
				else continue;
				proj_states->push_back(final_state);
			}

			iter_down <<= 1;
		}//END OF FOR L
		iter_up <<= 1;
	} //END OF FOR K
	//Remove duplicates
	std::sort(proj_states->begin(), proj_states->end());
	auto it = std::unique(proj_states->begin(), proj_states->end());
	proj_states->erase(it, proj_states->end());
}
sType create_min_k(int n_up, int n_down, hubbardParam* hubP);
int state_sym_locator(sType state, hubbardParam* hubP, std::vector<int>* k_index);
void Hepsilon(sType state, std::vector<sType>* proj_states, hubbardParam* hubP);
void u_jump_energy(sType right_state, Electrons elec, std::vector<sType>* states, std::vector<double>* energies, hubbardParam* hubP);
void epsilon_jump_energy(sType right_state, std::vector<sType>* states, std::vector<std::complex<double>>* energies, hubbardParam* hubP);
void epsilon_jump_energy(sType right_state, std::vector<sType>* states, std::vector<double>* energies, hubbardParam* hubP);
double state_energy(sType x, hubbardParam* hubP);
void calculate_epsilon_1d(hubbardParam* hubP);
void calculate_epsilon_3d(hubbardParam* hubP);

template <class A> void HuN_subspace_condition_expanding(A* sArr, uint64_t start, uint64_t end){
	/*******************************************************
	* Expends a subspace by applying the Hu hopping opoerator
	*
	* Parameters
	* ----------
	* sArr	: (A*) ptr to the StatesArr object to expand
	* start	: (uint64_t) first index to apply Ht on
	* end		: (uint64_t) index to stop apply Ht on
	*
	* Templates
	* ---------
	* A		: Any StatesArr child object
	*
	* Returns
	* -------
	* NONE
	********************************************************/
	for (uint64_t i = start; i < end; i++) {
		//HuN
		std::vector<sType> proj;
		HuN(sArr->get_at(i), &proj, sArr->sys_hubP.n_sites);
		for (uint64_t j = 0; j < proj.size(); j++) {
			sArr->add(proj.at(j));
		}
	}
}

template<class T, class U> void count_contribution_wH(
    T* fund, U* states, U* states_to_probe, double wH) {
	/*******************************************************
	* Sort the states according to their fund weight, creates a reduced
    * subspace keeping only the most dominant elements and write everything in
    * a .txt file
	*
	* Parameters
	* ----------
	* fund	: (std::vector<T>*) ptr to the fund vector
	* states: (U*) ptr to the StatesArr object to expand
	* sites	: (unsigned int) nbr of sites of the system
	* keep	: (double) weight to keep
	*
	* Templates
	* ---------
	* T		: double, std::complex<double>
	* U		: Any StatesArr child object
	*
	* Returns
	* -------
	* NONE
	********************************************************/
    sType size = states->get_length();
    wH *= pow(nrm2_blas(&size, fund),2);
	// Create index vector: [0, 1, 2, 3]
    std::vector<sType> indices(states->get_length());
    std::iota(indices.begin(), indices.end(), 0);

    // Sort indices based on weight of fundamental state
    std::sort(indices.begin(), indices.end(), [&](size_t i, size_t j) {
        return abs(fund[i]) > abs(fund[j]);
    });

	// Apply the sorting to weights and states
    //    sorted_fund[i] = fund[indices[i]];
    //    sorted_states[i] = states->get_at(indices[i]);
    //}

	//Removing old arrays
	//states->remove_all();
    //std::fill(fund, fund + states->get_length(), 0);

	//Creates a string to store the sorted fund
	std::string fund_txt = "";
	double cummul = 0;
	for (uLong i = 0 ; i < indices.size(); i++) {
        vType fund_index = fund[indices[i]];
        if constexpr(std::is_same_v<vType,double>) {
            cummul += (fund_index*fund_index);
        }
        else if constexpr(std::is_same_v<vType, std::complex<double>>) {
            cummul += (fund_index*std::conj(fund_index)).real();
        }

        //Add the states back in the States array but in order of their weight
        if (verbose > 9) std::cout<<cummul<<"/"<<wH<<std::endl;
        if (cummul < wH) {
            states_to_probe->add(states->get_at(indices[i]));
        }
		//Count how many states are needed to get to wH
		if (abs(cummul-wH) < 10e-8) {
            break;
			if(verbose > 99) {
                std::cout << cummul << "\tADD:"<<
                states->get_at(indices[i]) << "\t" << fund[indices[i]] << std::endl;
            }
		}

	}
}

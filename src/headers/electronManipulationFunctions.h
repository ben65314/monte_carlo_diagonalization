#pragma once
#include "Structures.h"


bool c_operator(sType* num, int index);
bool c_dag_operator(sType* num, int index);

void Ht(sType state, std::vector<sType>* proj_states, hubbardParam* hubP);

int Hu(sType state, unsigned char sites);

Electrons find_number_of_electron(sType state, unsigned char sites);

sType create_anti_ferro(unsigned int sites, int n_up, int n_down);

Electrons transform_NSz(int nElec, int spin);

void t_jump_energy(sType right_state, std::vector<sType>* states,
                   std::vector<double>* energies, hubbardParam* hubP);

double compute_mu(float mu, Electrons elec);

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
void u_jump_energy(sType right_state, Electrons elec, std::vector<sType>* states, std::vector<double>* energies, hubbardParam* hubP);
void epsilon_jump_energy(sType right_state, std::vector<sType>* states, std::vector<std::complex<double>>* energies, hubbardParam* hubP);
double state_energy(sType x, hubbardParam* hubP);
void calculate_epsilon_1d(hubbardParam* hubP);


template <class A> void Ht_subspace_condition_expanding(
        A* sArr, uint64_t start, uint64_t end){
	/*******************************************************
	* Expends a subspace by applying the Ht hopping opoerator
	*
	* Parameters
	* ----------
	* sArr	: (A*) ptr to the StatesArr object to expand
	* start	: (uint64_t) first index to apply Ht on
	* end	: (uint64_t) index to stop apply Ht on
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
		//Ht
		std::vector<sType> proj_Ht;
		Ht(sArr->get_at(i), &proj_Ht, &sArr->sys_hubP);
		//Add to sArr
		for (uint64_t j = 0; j < proj_Ht.size(); j++) {
			//Possible accept condition
			if (!accept_function(proj_Ht.at(j))) continue;
			sArr->add(proj_Ht.at(j));
		}
	}
}


template<class T, class U> void write_state_with_double(
    std::vector<T>* fund, U* states, unsigned int sites, double keep=1) {
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

	// Create index vector: [0, 1, 2, 3]
    std::vector<uLong> indices(fund->size());
    std::iota(indices.begin(), indices.end(), 0);

    // Sort indices based on keys
    std::sort(indices.begin(), indices.end(), [&](size_t i, size_t j) {
        return abs(fund->at(i)) > abs(fund->at(j));
    });

	// Apply the sorting to keys and values
    std::vector<T> sorted_fund(fund->size());
    std::vector<uLong> sorted_states(states->get_length());
    for (uLong i = 0; i < indices.size(); ++i) {
        sorted_fund[i] = fund->at(indices[i]);
        sorted_states[i] = states->get_at(indices[i]);
    }

	//Removing old arrays
	states->remove_all();
	fund->clear();

	//Creates a string to store the sorted fund
	std::string fund_txt = "";
	double cummul = 0;
	for (uLong i = 0 ; i < sorted_fund.size(); i++) {
		cummul += (double)(conjugate(sorted_fund.at(i))*sorted_fund.at(i)).real();
		fund_txt += to_string_pq(sorted_fund.at(i), 4, 14) + "\t"
            + to_string_pq((double)sorted_states.at(i), 10, 0) + "\t"
			+ to_string_pq((double)Hu(sorted_states.at(i), sites), 10, 0)+ "\t"
            + to_string_pq(cummul, 2, 16)+"\n";

		//Create reduced sampling size
		if (cummul < keep || keep==1) {
			if(verbose > 99) {
                std::cout << cummul << "\tADD:"<<
                sorted_states.at(i) << "\t" << sorted_fund.at(i) << std::endl;
            }
			states->add(sorted_states.at(i));
			fund->push_back(sorted_fund.at(i));
		}

	}

	//Writes
	char cwd[PATH_MAX];
	if (!getcwd(cwd, sizeof(cwd))) {
        std::cout << "Problem occured getting CWD" << std::endl;
    }
	std::string txtName = "/fund.txt";
	const std::string outFileName = cwd + txtName;
	std::ofstream outFile;
	//Write Q-matrix and ev
	outFile.open(outFileName);
	outFile<<fund_txt;
	outFile.close();


}

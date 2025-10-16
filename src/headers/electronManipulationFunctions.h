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

double compute_mu(float mu, Electrons elec);

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
		cummul += sorted_fund.at(i)*sorted_fund.at(i);
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

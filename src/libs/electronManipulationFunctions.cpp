#include "electronManipulationFunctions.h"
//Global variables
int verbose = 0;
//samplingParam global_sP;


bool c_operator(sType* num, int index){
	/*********************************************************
	* Construction operator of the Fock formalism
	*
	* Parameters
	* ----------
	* num		: (sType*) state number to change
	* index		: (int) where to create
	*
	* Returns
	* -------
	* c_success	: true if the creation was successful
	***********************************************************/
	sType cd_op = 1UL << index;
	bool c_success = true;
	//if the bit operation is the same as the addition the creation is good
	if((*num | cd_op) != *num){
		*num += cd_op;
	}
	else{
		c_success =  false;	
	}
	return c_success;
}
bool cDag_operator(sType* num, int index){
	/*********************************************************
	* Destruction operator of the Fock formalism
	*
	* Parameters
	* ----------
	* num		: (sType*) state number to change
	* index		: (int) where to create
	*
	* Returns
	* -------
	* d_success	: true if the destruction was successful
	***********************************************************/
	sType cd_op = 1UL << index;
	bool d_success = true;
	//if the bit operation is the same as the substraction the creation is good
	if ((*num & cd_op) != 0){
		*num -= cd_op;
	}
	else{
		d_success = false;	
	}
	return d_success;
}
void Ht(sType state, std::vector<sType>* proj_states, hubbardParam* hubP) {
	/****************************************************************
	* Applies the hopping operator of the Hubbard Hamiltonian
	*
	* Parameters
	* ----------
	* state			: (sType) initial state in the Fock formalism
	* proj_states	: (std::vector<sType>*) Array of the possible accessible states after the Hamiltonian
	* hubP			: (hubbardParam*) System parameters
	*
	* Templates:
	* ----------
	* T				: int, short, long, unsigned
	*
	* Returns
	* -------
	* NONE
	******************************************************************/
	//States after applied Hamiltonian
	unsigned char sites = hubP->n_sites;
	sType one = 1;
	proj_states->reserve(2 * sites);

	const sType state_num = state;
	sType from_down = one << (sites - 1), from_up = one << (2 * sites - 1);
	for(unsigned char i = 0; i < sites; i++){
		sType to_down = one << (sites - 1), to_up = one << (2 * sites - 1);
		for(unsigned char j = 0; j < sites; j++){
			if (hubP->tMatrix[i * sites + j] == 0){
				to_down >>= 1;
				to_up >>= 1;
				continue;
			}
			
			if (((from_down & state_num) != 0) && ((to_down & state_num ) == 0)) {
				proj_states->push_back((state_num | to_down) ^ from_down);
			}
			
			if (((from_up & state_num) != 0) && ((to_up & state_num ) == 0)) {
				proj_states->push_back((state_num | to_up) ^ from_up);
			}
			
			to_down >>= 1;
			to_up >>= 1;
		}
		from_down >>= 1;
		from_up >>= 1;
	}
}

int Hu(sType state, unsigned char sites) {
	/**********************************************************
	* Computes the number of coulomb interaction present in the state
	*
	* Parameters
	* ----------
	* state : (sType) state to check interatction
	* sites : (unsigned char) number of sites of the system
	*
	* Returns
	* -------
	* nU	: (int) number of sites that have a double occupation
	***********************************************************/
	sType same_occupation = state & (state << sites); 
	int nU = oneCounter(same_occupation);
	return nU;
}
Electrons findNumberOfElectron(sType state, unsigned char sites) {
	/*****************************************************************
	* Counts the number of electron of each spin for a given state
	*
	* Parameters
	* ----------
	* state : (sType) state to count electrons
	* sites : (unsigned char) number of sites of the state
	*
	* Returns
	* -------
	* elec	: (Electrons) countains the up and down electrons
	******************************************************************/
	sType up_counter = 1UL << sites;
	sType down_counter = 1;

	Electrons elec;
	
	for (unsigned char j = 0; j < sites; j++){
		elec.up += ((state & up_counter) != 0);
		elec.down += ((state & down_counter) != 0);

		up_counter <<= 1;
		down_counter <<= 1;
	}

	return elec;
}

sType createAntiFerro(unsigned int sites, int n_up, int n_down){
	/*******************************************************************
	* Creates the best antiFerro state
	*
	* Parameters
	* ----------
	* sites : (unsigned int) Number of sites of the system
	* n_up	: (int) number of up electrons
	* n_down : (int) number of down electrons
	*
	* Returns
	* -------
	* bi_state : (sType) Anti-Ferro state 
	*****************************************************************/
	sType bi_state = 0;
	unsigned long adder = 1;
	unsigned long one = 1;
	//Creates antiFerro State if the number of state is even
	for(int i = 0; i < n_down; i++){
		if(adder >= (sType)(one << sites))
		{
			adder = 2;
		}
		bi_state |= adder;
		adder <<= 2;
	}

	adder = one << (sites + 1);
	for(int i = 0; i < n_up; i++){
		if(adder >= (one << (2 * sites)))
		{
			adder = one << sites;
		}
		bi_state |= adder;
		adder <<= 2;
	}

	return bi_state;
}

//void createEachFrequency(hubbardParam* hubP, int nUp, int nDown, std::vector<T>* states)

Electrons transformNSz(int n_elec, int spin) {
	/*****************************************
	* Transforms the total number of electron of the system and the total spin in an Electrons structure countaining the number of ups and downs
	*
	* Parameters
	* ----------
	* n_elec: (int) Total number of electrons of the system
	* spin	: (int) Total spin of the system
	* 
	* Returns
	* -------
	* elec	: (Electrons) countains the number of spin up and down
	****************************************/
	Electrons elec;
	elec.up = (spin + n_elec)/2;
	elec.down = (n_elec - spin)/2;
	return elec;
}

//Jump ENERGIES
void tJumpEnergy(sType right_state, std::vector<sType>* states, std::vector<double>* energies, hubbardParam* hubP) {
	/*******************************************************************
	* Calculates the energy of a t jump between two given states
	*
	* Parameters
	* ----------
	* right_state	: (sType) state to jump from
	* states		: (std::vector<sType>*) receptacles of the states accessible from right_state
	* energies		: (std::vector<double>*) receptacles of the energiesfor each states
	* hubP		: (hubbardParam*) System parameters
	*
	* Returns
	* -------
	* NONE
	*****************************************************************/
	//Calculates the electron jump sites
	int sites = hubP->n_sites;
	for (int i = 0; i < sites; i++) {//From
		for (int j = i+1; j < sites; j++) {//To
			//Jump energy
			double jumpFactor = hubP->tMatrix.at(i * sites + j);
			if (jumpFactor == 0) continue; 
			
			for (int k = 0; k < 2; k++) {//Iteration over spins
				sType temp_stateJI = right_state;
				sType temp_stateIJ = right_state;

				int indexJ = (sites-j-1 + k*sites);
				int indexI = (sites-i-1 + k*sites);
				//From J -> I
				if (cDag_operator(&temp_stateIJ,indexJ) && c_operator(&temp_stateIJ,indexI)) {
					//Add to the receptacle
					states->push_back(temp_stateIJ);

					//Phase
					temp_stateIJ >>= indexJ + 1;
					char phase = 1;
					for (int l = indexJ + 1; l < indexI; l++) {
						if ((temp_stateIJ & 1) == 1) phase *= -1;  
						temp_stateIJ >>= 1;	
					}
					//Add the energy to the receptacle
					energies->push_back(jumpFactor * phase);
				}	
				//From I -> J
				if (cDag_operator(&temp_stateJI,indexI) && c_operator(&temp_stateJI,indexJ)) {
					//Add to the receptacle
					states->push_back(temp_stateJI);

					//Phase
					temp_stateJI >>= indexJ + 1;
					char phase = 1;
					for (int l = indexJ + 1; l < indexI; l++) {
						if ((temp_stateJI & 1) == 1) phase *= -1;  
						temp_stateJI >>= 1;	
					}
					//Add the energy to the receptacle
					energies->push_back(jumpFactor * phase);
				}
			}
		}
	}
}

//void Ht_subspace_condition_expanding(A* sArr, uint64_t start, uint64_t end)
//void writeStateWithDouble(std::vector<T>* fund, U* states, unsigned int sites, double keep=1)

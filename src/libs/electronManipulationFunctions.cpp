#include "electronManipulationFunctions.h"
#include "Structures.h"
#include "blas_lapack_wrappers.h"
#include "basicFunctions.h"


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
	if ((*num | cd_op) != *num) {
		*num += cd_op;
	}
	else {
		c_success =  false;
	}
	return c_success;
}
bool c_dag_operator(sType* num, int index){
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
	if ((*num & cd_op) != 0) {
		*num -= cd_op;
	}
	else {
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
	* proj_states	: (std::vector<sType>*) Array of the possible accessible
    *                                       states after the Hamiltonian
	* hubP			: (hubbardParam*) System parameters
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
    //Where the electron is
	sType from_down = one << (sites - 1);
    sType from_up = one << (2 * sites - 1);

	for(unsigned char i = 0; i < sites; i++){
        //Where the electron is going
		sType to_down = one << (sites - 1);
        sType to_up = one << (2 * sites - 1);

		for(unsigned char j = 0; j < sites; j++){
            //Skips the iteration if the jump has no amplitude in the t matrix
			if (hubP->t_matrix[i * sites + j] == 0){
				to_down >>= 1;
				to_up >>= 1;
				continue;
			}

            //Down electron jump
			if (((from_down & state_num) != 0)
                && ((to_down & state_num ) == 0)) {
				proj_states->push_back((state_num | to_down) ^ from_down);
			}

            //Up electron jump
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
	int nU = one_counter(same_occupation);
	return nU;
}
Electrons find_number_of_electron(sType state, unsigned char sites) {
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
		elec.up += (state & up_counter) != 0;
		elec.down += (state & down_counter) != 0;

		up_counter <<= 1;
		down_counter <<= 1;
	}

	return elec;
}

sType create_anti_ferro(unsigned int sites, int n_up, int n_down){
	/*******************************************************************
	* Creates the best antiferromagnetic state
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

Electrons transform_NSz(int n_elec, int spin) {
	/*****************************************
	* Transforms the total number of electron of the system and the total spin
    * in an Electrons structure countaining the number of ups and downs
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
void t_jump_energy(sType right_state, std::vector<sType>* states,
                   std::vector<double>* energies, hubbardParam* hubP) {
	/*******************************************************************
	* Calculates the energy of a t jump between two given states
	*
	* Parameters
	* ----------
	* right_state	: (sType) state to jump from
	* states		: (std::vector<sType>*) receptacles of the states
    *                                       accessible from right_state
	* energies		: (std::vector<double>*) receptacles of the energiesfor
    *                                        each states
	* hubP		    : (hubbardParam*) System parameters
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
			double jumpFactor = hubP->t_matrix.at(i * sites + j);
			if (jumpFactor == 0) continue;

			for (int k = 0; k < 2; k++) {//Iteration over spins
				sType temp_stateJI = right_state;
				sType temp_stateIJ = right_state;

				int indexJ = (sites-j-1 + k*sites);
				int indexI = (sites-i-1 + k*sites);
				//From J -> I
				if (c_dag_operator(&temp_stateIJ, indexJ)
                    && c_operator(&temp_stateIJ, indexI)) {
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
				if (c_dag_operator(&temp_stateJI, indexI)
                    && c_operator(&temp_stateJI, indexJ)) {
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

//void writeStateWithDouble(std::vector<T>* fund, U* states, unsigned int sites, double keep=1)

double compute_mu(float mu, Electrons elec){
	/**************************************************
	* Computes mu interaction used in matrixCreation
	*
	* Parameters
	* ----------
	* mu: (float) mu value
	* elec: (Electrons) number of electrons of the state
	*
	* Returns
	* ------ -
	* result: (double) mu value
	************************************************/

	double result;
	//Counts number of electrons
	long n_electron = elec.up + elec.down;
	//Computs chemical potential
	result = -n_electron * mu;
	return result;
}

// k-basis functions
//

void u_jump_energy(sType right_state, Electrons elec, std::vector<sType>* states, std::vector<double>* energies, hubbardParam* hubP) {
	/*******************************************************************
	* Calculates the energy of a U_N jump between two given states
	*
	* Parameters
	* ----------
	* right_state		: (long) state to jump from
	* elec				: (Electrons) number of electrons of the right_state
	* states_accessible : (std::vector<sType>) states to project to
	* energies			: (std::vector<double>) energies for states to project to
	* hubP			: (hubbardParam*) System parameters
	*
	* Returns
	* -------
	* NONE
	*****************************************************************/
	double jump_energy = hubP->u / hubP->n_sites;
	double stay_energy = jump_energy * elec.up * elec.down;

	//No mouvement
	states->push_back(right_state);
	energies->push_back(stay_energy);

	sType iter_up = 1UL << hubP->n_sites;
	//Checks for all k
	for (int k = hubP->n_sites - 1; k >= 0; k--) {
		//Check for all destruction sites
		sType iter_down = 1L;
		for (int l = hubP->n_sites - 1; l >= 0; l--) {
			sType evolving_state = right_state;
			//long change = 0;

			if (((evolving_state & iter_up) != 0) && ((evolving_state & iter_down) != 0)) {

				evolving_state ^= (iter_up | iter_down);
				//change -= (iterUp + iterDown);
			}
			else {iter_down <<= 1;	continue;}
			for (int q = 1; q < hubP->n_sites; q++) {
				long final_state;
				int qk = q+k;
				long where_up;
				if (qk >= hubP->n_sites) {
					where_up = iter_up << (hubP->n_sites - q);
					qk -= hubP->n_sites;
				}
				else where_up = iter_up >> (q);

				int ql = l - q;
				long where_down;
				if (ql < 0) {
					where_down = iter_down >> (hubP->n_sites - q);
					ql += hubP->n_sites;
				}
				else where_down = iter_down << q;


				if (((evolving_state | where_up) != evolving_state) && ((evolving_state | where_down) != evolving_state)) {
					final_state = evolving_state | (where_up | where_down);
					//change += (whereUp + whereDown);
				}
				else continue;

				states->push_back(final_state);
				char phase = 1;
				int start_down, start_up;
				if (qk < k) start_up = k;
				else {
					start_up = qk;
					qk = k;
				}

				if (ql < l) start_down = l;
				else {
					start_down = ql;
					ql = l;
				}
				sType down_state = final_state >> (hubP->n_sites - start_down);

				for (int m = start_down - 1; m > ql; m--) {
					if ((down_state & 1) == 1) phase *= -1;
					down_state >>= 1;
				}

				sType up_state = final_state >> (2 * hubP->n_sites - start_up);
				for (int m = start_up - 1; m > qk; m--) {
					if ((up_state & 1) == 1) phase *= -1;
					up_state >>= 1;
				}
				energies->push_back(jump_energy * phase);
			}//END OF FOR q

			iter_down <<= 1;
		}//END OF FOR L
		iter_up <<= 1;
	} //END OF FOR K
}


void epsilon_jump_energy(sType right_state, std::vector<sType>* states, std::vector<std::complex<double>>* energies, hubbardParam* hubP) {
	/******************************************************
	* Calculates the energy of a He jump between two given states
	*
	* Parameters
	* ----------
	* right_state	: (sType) state to jump from
	* states		: (std::vector<sType>*) receptacles of the states accessible from right_state
	* energies		: (std::vector<std::complec<double>>*) receptacles of the energies for each states
	* hubP		: (hubbardParam*) System parameters
	*
	* Returns
	* -------
	* NONE
	*******************************************************/
	//States after applied Hamiltonian
	int sites = hubP->n_sites;
	const sType state_num = right_state;
	states->push_back(right_state);
	energies->push_back(0);

	sType from_down = 1UL << (sites - 1), from_up = 1UL << (2 * sites - 1);

	for(int i = 0; i < sites; i++){//From
		sType to_down = 1UL << (sites - 1), to_up = 1UL << (2 * sites - 1);

		for(int j = 0; j < sites; j++){//To
			if (hubP->matEpsilon[i * sites + j]==std::complex<double>(0,0)){
				to_down >>= 1;
				to_up >>= 1;
				continue;
			}
			//Is there an electron there and can we jump there
			bool has_down, has_up, can_up, can_down, standby_up, standby_down;
			has_down = ((from_down & state_num) != 0);
			can_down = ((to_down & state_num ) == 0);

			has_up = ((from_up & state_num) != 0);
			can_up = ((to_up & state_num ) == 0);

			standby_up = (from_up == to_up);
			standby_down = (from_down == to_down);

			//Spin down
			if ((has_down && can_down) || (has_down && standby_down)) {
				sType new_state = (state_num+to_down-from_down);

				if (new_state == right_state) {
					(*energies)[0] += hubP->matEpsilon[i * sites + j];
				}
				else {
					states->push_back(new_state);

					char phase = 1;
					char start = (i > j) ? i : j;
					char end = (i > j) ? j : i;
					new_state >>= (hubP->n_sites - start);
					for (int l = start - 1; l > end; l--) {
						if ((new_state & 1) == 1) phase *= -1;
						new_state >>= 1;
					}
					energies->push_back(hubP->matEpsilon[i * sites + j] * (double)phase);
				}

			}
			//Spin up
			if ((has_up && can_up) || (has_up && standby_up)) {
				sType new_state = (state_num+to_up-from_up);

				if (new_state == right_state) {
					(*energies)[0] += hubP->matEpsilon[i * sites + j];
				}
				else {
					states->push_back(new_state);

					char phase = 1;
					char start = (i > j) ? i : j;
					char end = (i > j) ? j : i;
					new_state >>= (hubP->n_sites - start);
					for (int l = start - 1; l > end; l--) {
						if ((new_state & 1) == 1) phase *= -1;
						new_state >>= 1;
					}
					energies->push_back(hubP->matEpsilon[i * sites + j] * (double)phase);
				}
			}
			to_down >>= 1;
			to_up >>= 1;
		}
		from_down >>= 1;
		from_up >>= 1;
	}
}

double state_energy(sType x, hubbardParam* hubP){
	/*****************************************
	* Calculates the energy <x|H|x> of a given state
	*
	* Parameters
	* ----------
    * x     : (sType) state
	* hubP	: (hubbardParam) System parameters
	*
	* Returns
	* -------
	* energy: (double) energy of the state
	****************************************/

    double energy = 0;
    //Epsilon terms
    sType scan = 1 << 2*hubP->n_sites;
    for (int i=0; i < 2*hubP->n_sites; i++){
        scan >>= 1;
        if ((x & scan) != 0){
            energy += hubP->matEpsilon[(i%hubP->n_sites)*(1+hubP->n_sites)].real();
        }
    }

    //U terms
    Electrons elec = find_number_of_electron(x, hubP->n_sites);
	energy += hubP->u / hubP->n_sites * elec.up * elec.down;

    return energy;
}

void calculate_epsilon_1d(hubbardParam* hubP){
	/*****************************************
	* Calculates the epsilon values matrix in a 1D lattice
	*
	* Parameters
	* ----------
	* hubP	: (hubbardParam) System parameters
	*
	* Returns
	* -------
	* NONE
	****************************************/
	std::complex<double> value = 0;
	for (int i = 0; i < hubP->n_sites; i++) {
		for (int j = 0; j < hubP->n_sites; j++) {
			value = 0;
			for (int k = 0; k < hubP->n_sites; k++){
				for (int l = 0; l < hubP->n_sites; l++){
					double delta = 2 * M_PI * (double)(-(double)(i * k) + j * l);
					delta /= hubP->n_sites;
					value += hubP->t_matrix[k * hubP->n_sites + l] * exp(std::complex<double>(0,1) * delta) / (double) hubP->n_sites;
				}//END OF FOR L
			}//END OF FOR K
			hubP->matEpsilon[i * hubP->n_sites + j] = remove_zeros(value);

		}//END OF FOR J
	}//END OF FOR I
}


int INC = 1;
void calculate_epsilon_3d(hubbardParam* hubP){
	/*****************************************
	* Calculates the epsilon values matrix in a 1D lattice
	*
	* Parameters
	* ----------
	* hubP	: (hubbardParam) System parameters
	*
	* Returns
	* -------
	* NONE
	****************************************/
    sType* dim = &(hubP->DIM);
	std::complex<double> value = 0;
	for (int k = 0; k < hubP->n_sites; k++) {
		for (int q = 0; q < hubP->n_sites; q++) {
			value = 0;
			for (int i = 0; i < hubP->n_sites; i++){
				for (int j = 0; j < hubP->n_sites; j++){
                    double delta =ddot_(dim,hubP->K.data()+k*(*dim),
                                        &INC,hubP->R.data()+i*(*dim),&INC) - 
                                  ddot_(dim,hubP->K.data()+q*(*dim),
                                        &INC,hubP->R.data()+j*(*dim),&INC);

					value += hubP->t_matrix[i * hubP->n_sites + j] * exp(std::complex<double>(0,1) * delta) / (double) hubP->n_sites;
				}//END OF FOR J
			}//END OF FOR I
			hubP->matEpsilon[k * hubP->n_sites + q] = remove_zeros(value);

		}//END OF FOR Q
	}//END OF FOR K
}
void Hepsilon(sType state, std::vector<sType>* proj_states, hubbardParam* hubP) {
	/****************************************************************
	* Applies the hopping operator of the Hubbard Hamiltonian
	*
	* Parameters
	* ----------
	* state			: (sType) initial state in the Fock formalism
	* proj_states	: (std::vector<sType>*) Array of the possible accessible
    *                                       states after the Hamiltonian
	* hubP			: (hubbardParam*) System parameters
	*
	* Returns
	* -------
	* NONE
	******************************************************************/
	//States after applied Hamiltonian
	unsigned char sites = hubP->n_sites;
	sType one = 1;

	const sType state_num = state;
    //Where the electron is
	sType from_down = one << (sites - 1);
    sType from_up = one << (2 * sites - 1);

	for(unsigned char i = 0; i < sites; i++){
        //Where the electron is going
		sType to_down = one << (sites - 1);
        sType to_up = one << (2 * sites - 1);

		for(unsigned char j = 0; j < sites; j++){
            //Skips the iteration if the jump has no amplitude in the t matrix
			if (hubP->matEpsilon[i * sites + j] == (std::complex<double>)0){
				to_down >>= 1;
				to_up >>= 1;
				continue;
			}

            //Down electron jump
			if (((from_down & state_num) != 0)
                && ((to_down & state_num ) == 0)) {
				proj_states->push_back((state_num | to_down) ^ from_down);
			}

            //Up electron jump
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

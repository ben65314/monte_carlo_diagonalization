#include "basicFunctions.h"

struct hubbardParam;

void isSuccess(bool status) {
	/*****************************************
	* Print the success value with a change of colors, for benchmarks
	*
	* Parameters
	* ----------
	* status : (bool) true if a success wants to be print false if else
	*
	* Returns
	* -------
	* NONE
	*****************************************/
	std::string result;
	if (status) {
		result = "SUCCESS";
		std::cout << "\033[1;32m~SUCCESS\033[0m\n";
	}
	else {
		result = "FAIL";
		std::cout << "\033[1;6;31m~FAIL\033[0m\n";
	}	
}
void cct(std::string text, int color) {
    /*****************************************
	* Print the text in a specified color
	*
	* Parameters
	* ----------
	* text	: (string) text to change color
	* color : (int) color value
	*
	* Returns
	* -------
	* NONE
	*****************************************/
	std::cout << "\033[1;" << color << "m" << text << "\033[0m";	
}
int getNumberOfBlocks(int sites) {
	/*****************************************
	* Counts the number of different blocks the Hamiltonian will divide into
	*
	* Parameters
	* ----------
	* sites		: (int) Number of sites of the system
	*
	* Returns
	* -------
	* nbrBlocks : (int) Number of different blocks
	* ****************************************/
	int nbrBlocks = sites + 1;
	nbrBlocks *= nbrBlocks;
	return nbrBlocks;
}
double findMinOfArray(std::vector<double> array) {
	/*****************************************
	* Finds the minimum value of an array by just looking every element
	* 
	* Parameters
	* ----------
	* array : (std::vector<double>) array to search through
	*
	* Returns
	* -------
	* min	: (double) minimum value countained in the array
	*****************************************/
	double min = array.at(0);
	for (unsigned long i = 1; i < array.size(); i++) {
		if (min > array.at(i)) {
			min = array.at(i);
		}
	}
	return min;
}
unsigned long comb(uLong n, uLong r) {
	/*****************************************
	* Calculates the combination term nCr
	*
	* Parameters
	* ----------
	* n		: (unsigned long) n value of the formula
	* r		: (unsigned long) r value of the formula
	*
	* Returns
	* -------
	* result: (unsigned long) possible combinations 
	****************************************/
	unsigned long num = 1;
	unsigned long den = 1;
	if ((n - r) > r) {
		for (uLong i = n; i > (n - r); i--) num *= i;

		for (uLong i = r; i > 1; i--) den *= i;
	}
	else {
		for (uLong i = n; i > r; i--) num *= i;
		
		for (uLong i = (n-r); i > 1; i--) den *= i;
	}
	unsigned long result = num / den;
	if (n < 1 && r > 0) result = 0 ;
	return result;
}
unsigned long combSpecified(uLong nU, uLong sites, uLong up, uLong down) {
	/*****************************************
	* Calculates the combination term for a specified number of up electrons, down electrons and double occupation
	*
	* Parameters
	* ----------
	* nU	: (unsigned long) number of double occupation
	* sites : (unsigned long) number of sites
	* up	: (unsigned long) number of ups
	* down	: (unsigned long) number of downs
	*
	* Returns
	* -------
	* result: (unsigned long) possible combinations 
	****************************************/
	if (up < nU || down < nU) {
		std::cout<<"WARNING - number of electrons inferior to number of double occupation"<<std::endl;
		return 0;
	}

	up -= nU;
	down -= nU;
	uLong result = comb(sites,nU) * comb(sites-nU,up) * comb(sites-nU-up,down);

	return result;
}
double fermiDiracFunction(float x, float beta, float mu) {
	/*****************************************
	* Evaluates the point x of a a Fermi-Dirac distribution with parameters beta and mu given
	*
	* Parameters
	* ----------
	* x		: (float) point to evaluate the Fermi-Dirac distribution
	* beta	: (float) beta coefficient in the distribution
	* mu	: (float) mu coefficient in the distribution (x=mu -> 0.5)
	*
	* Returns
	* -------
	* n		: (double) value of the Fermi-Dirac distribution at x
	******************************************/
	double exp_ =  exp(beta*(x-mu));
	double n = 1/(exp_ + 1);

	return n;
}
double boltzmannDistributionFunction(float dE, float beta) {
	/*****************************************
	* Evaluates the probability of change of the energy dE in a Boltzmann distribution with parameter beta given
	*
	* Parameters
	* ----------
	* dE	: (float) point to evaluate the Boltzmann distribution
	* beta	: (float) beta coefficient in the distribution
	*
	* Returns
	* -------
	* exp_	: (double) value of the Boltzmann distribution at dE
	****************************************/
	double exp_ =  exp(-beta*dE);

	return exp_;
}

double calculateSD(std::vector<double> data) {
	/*****************************************
	* Calculates the standard deviation of a data sample
	*
	* Parameters
	* ----------
	* data	: (std::vector<double>) data sampled to compute
	*
	* Returns
	* -------
	* std	: (double) standard deviation of the sample
	****************************************/
	double sum = 0.0, mean, standardDeviation = 0.0;
	int i;
	int len = data.size();
	for (i = 0; i < len; ++i) {
		sum += data.at(i);
	}

	mean = sum / len;

	for (i = 0; i < len; ++i) {
		double temp = data.at(i) - mean;
		standardDeviation += temp*temp;
	}

	double std = sqrt(standardDeviation / len);
	return std;
}
int degFundamentalCheck(double* eV, uLong n, double eps) {
	/*********************************************************
	* Checks the degeneracy level of the fundamental
	*
	* Parameters
	* ----------
	* eV	: (double*) eigen values computed 
	* n		: (unsigned long) number of ev
	* eps	: (double) precision for degeneracy 
	*
	* Returns
	* -------
	* NONE
	*******************************************************/
	int deg_degree = 1;
	for (uLong i = 0; i < n; i++) {
		//Increase counter if energies are the same
		if (abs(eV[i]-eV[i+1]) < eps) deg_degree++;
		else break;
	}
	return deg_degree;
}

//void makeIdentity(T* matrix, R rows) in header file
//void makeTridiag(T* matrix, R* diag, R* h_diag, U size) in header file
//void normalize(double* vec, R size) in header file
//void initialVector(R const SIZE, T* v, uInt seed) in header file
//void gramSchmidtIteration(double* prev_vectors, R const size_vectors, R const num_vectors, double* to_ortho_vector) in header file

std::string toStringP(const double a_value, const int n) {
	/*****************************************
	* Transforms a double value in a string with a given precicion
	*
	* Parameters
	* ----------
	* a_value	: (const double) value to stringify
	* n			: (const int) precision requested
	*
	* Returns
	* -------
	* string	: (double) converted to string with precision n
	***************************************/
	std::ostringstream out;
	out.precision(n);
	out << std::fixed << a_value;
	return std::move(out).str();
}

template <> std::string toStringP_Q(const double a_value, const uInt integers, const uInt decimals) {
	/*****************************************
	* Transforms a double value in a string with a given precicion
	*
	* Parameters
	* ----------
	* a_value	: (const double) value to stringify
	* integers	: (const unsigned int) precision of integer part
	* decimals	: (const unsigned int) precision of decimal part
	*
	* Returns
	* -------
	* string	: (double) converted to string
	***************************************/
	char s[100];
	if (abs(a_value) > 10e-10) {
		if (a_value < 0) sprintf(s,"% *.*f",integers+decimals,decimals,a_value);
		else sprintf(s,"% *.*f",integers+decimals,decimals,a_value);
	}
	else sprintf(s,"% *.*f",integers+decimals,decimals,0.0);

	std::string str(s);
	return str;
}

template <> double removeZeros(double a){
	/*****************************************
	if a value of a real number is near zero, will make it zero

	Parameters
	----------
	a: (double) number to verify

	Returns
	-------
	a: (double) number with near zero ->  zero
	****************************************/
	if (abs(a) < 10e-13){a = 0;}	
	return a;
}

//std::string writeVector(const double* vec, R size) in header file
//std::string writeMatrix(const double* mat, R rows, R cols, int integers, int precision) in header file
//std::string printVector(const T* vec, R size) in header file
//std::string printMatrix(const T* mat, R rows, R cols, int integers, int precision) in header file

std::string timeFormating(std::chrono::time_point<std::chrono::high_resolution_clock> start, std::chrono::time_point<std::chrono::high_resolution_clock> end){
	/*****************************************
	* Takes two time point and puts them in a string accordingly to their size
	*
	* Parameters
	* ----------
	* start	: (time_point) First time point
	* end	: (time_point) Second time point

	* Returns
	* -------
	* time	: (string) time passed between the two points in string format
	****************************************/
	std::string time="";
	double min, sec, ms, mus;
	min = sec = ms = mus = 0;

	//Minutes
	if(std::chrono::duration_cast<std::chrono::minutes>(end - start).count() > 1){
		min = floor((double)std::chrono::duration_cast<std::chrono::minutes>(end - start).count());
		time += toStringP(min, 0) + "min ";
	}
	//Seconds
	if(std::chrono::duration_cast<std::chrono::seconds>(end - start).count() > 1){
		sec = floor((double)std::chrono::duration_cast<std::chrono::seconds>(end - start).count());
		time += toStringP(sec- min * 60, 0) + "s ";
	}
	//Milliseconds
	if(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() > 1){
		ms = floor((double)std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
		time += toStringP(ms - sec * 1000, 0) + "ms ";
	}
	//Microseconds
	if(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() > 1){
		mus = floor((double)std::chrono::duration_cast<std::chrono::microseconds>(end - start).count());
		time += toStringP(mus - 1000 * ms, 0) + "μs ";
	}
	return time;
}

unsigned int oneCounter(sType num) {
	/*****************************************
	* Counts the number of ones in a given number in a binary representation
	*
	* Parameters
	* ----------
	* num	: (sType) pointer of the complex matrix to stringify
	* 
	* Returns
	* -------
	* count : (unsigned int) number of ones.
	****************************************/
	sType one = 1;
	unsigned int n_elements = sizeof(sType) * 8;
	unsigned int count = 0;
	for (unsigned int i = 0; i < n_elements; i++) {
		if ((one & num) != 0) {count++;}
		one <<= 1;
	}

	return count;
}
sType rotr(sType bitVec, unsigned char rotIndex, unsigned char size) {
	/*****************************************
	* Shifts bits to the right
	*
	* Parameters
	* ----------
	* bitVec	: (sType) number to rotr
	* rotIndex	: (unsigned char) move to right number of times
	* size		: (unsigned char) number of bit string
	* 
	* Returns
	* -------
	****************************************/
	return (((bitVec >> rotIndex) | (bitVec << (size - rotIndex))) & ((sType)pow(2,size)-1));
}

void combinationRecursive(std::vector<uShort> empty_spaces, sType current_num, uShort left_to_place, sType* all_comb, sType* placed, sType* it){
	/****************************************************************
	* Creates the combinations of the number specified <leftToPlace> to each place accessible <empty_spaces>. This function is specially used when creating double occupation states.
	* 
	* Parameter
	* ---------
	* empty_spaces	: (std::vector) Emplacement where we can place an electron
	* current_num	: (unsigned long) coutains the previously placed electron
	* left_to_place	: (unsigned short) number of electrons to placed
	* all_comb		: (unsigned long*) Countains all the possible spaces
	* placed		: (unsigned long*) number of state placed
	*
	*
	* Returns
	* -------
	* NONE
	*****************************************************************/
	
	if (left_to_place == 0){
		(*it) += 1;
		*(all_comb + *placed) = current_num;
		(*placed) += 1;
		return;
	}
	
	for (uInt i = 0; i < empty_spaces.size(); i++) {
		uLong added_num = current_num ^ (1ul << empty_spaces.at(i));
		std::vector<uShort> updated_empty_spaces = empty_spaces;
		updated_empty_spaces.erase(updated_empty_spaces.begin(), updated_empty_spaces.begin() + i + 1);


		combinationRecursive(updated_empty_spaces, added_num, left_to_place - 1, all_comb, placed, it);
	}
}

void combinationDoubleOccupation(uShort N, uShort sites, std::vector<sType>* all_double_occupation) {
	/****************************************************************
	* Creates all the combinations of states with <N> double electron occupation in a system of <sites> sites. This will coutain only the <N> double occupation without any single electron occupation on any site.
	* 
	* Parameter
	* ----------
	* N						: (uShort) Number of required double occupation
	* sites					: (uShort) Number of sites of the system
	* allDoubleOccupation	: (std::vector) all the <N> double occupation states
	*
	* Returns
	* -------
	* NONE
	*****************************************************************/
	uint64_t combU = comb(sites,N);
	*all_double_occupation = std::vector<sType>(combU,0);
	uLong placed = 0;

	//All sites accessible
	std::vector<uint16_t> chosenSites(sites);
	std::iota(chosenSites.begin(), chosenSites.end(),0);
	uint64_t current_num = 0;

	//Recursive combination function
	//StatesArr* combinations = new StatesArr(combU);
	uLong it = 0;
	combinationRecursive(chosenSites,current_num,N,all_double_occupation->data(), &placed, &it);

	//Shift the states found to create the double occupation
	for (unsigned int i = 0; i < all_double_occupation->size(); i++) {
		sType temp = all_double_occupation->at(i);
		all_double_occupation->at(i) = ((temp << sites) | temp);
	}
}

void combinationRecursiveAddingSingle(std::vector<uShort> empty_spaces_up, std::vector<uShort> empty_spaces_down, sType current_num, uShort n_up, uShort n_down, uShort sites, sType* all_comb, sType* placed) {
	/****************************************************************
	* Completes the <currentNum> with <nUp> and <nDown> single occupationelectrons
	* 
	* Parameter
	* ----------
	* empty_spaces_up	: (std::vector) Emplacement where we can place an up electron
	* empty_spaces_down	: (std::vector) Emplacement where we can place a down electron
	* current_num		: (ulong) coutains the previously placed electrons
	* n_up				: (uShort) number of up electrons to place
	* n_down			: (uShort) number of down electrons to place
	* sites				: (uShort) number of sites of the system
	* all_comb			: (sType*) Countains all the possible spaces
	* placed			: (unsigned long*) Number of elements placed in allComb
	*
	* Returns
	* -------
	* NONE
	*****************************************************************/
	//If all the electrons are placed.
	if (n_up <= 0 && n_down <= 0){
		*(all_comb + *placed) = current_num;
		(*placed) += 1;
		return;
	}
	uLong added_num = current_num;

	if (n_up > 0) {	//Add an up electron
		for (uint i = 0; i < empty_spaces_up.size(); i++) {
			added_num = current_num ^ (1ul << (empty_spaces_up.at(i) + sites));

			std::vector<uShort> updated_empty_spaces_up = empty_spaces_up;
			std::vector<uShort> updated_empty_spaces_down = empty_spaces_down;
			//Remove the access for down electrons to be put here
			updated_empty_spaces_down.erase(std::remove(updated_empty_spaces_down.begin(), updated_empty_spaces_down.end(), empty_spaces_up.at(i)), updated_empty_spaces_down.end());
			//Remove the access for up electrons to be put here
			updated_empty_spaces_up.erase(updated_empty_spaces_up.begin(), updated_empty_spaces_up.begin() + i + 1);

			//Recursive
			combinationRecursiveAddingSingle(updated_empty_spaces_up,updated_empty_spaces_down, added_num, n_up-1, n_down, sites, all_comb, placed);
		}

	}
	else if (n_down > 0) {	//Add a down electron
		for (uint i = 0; i < empty_spaces_down.size(); i++) {
			added_num = current_num ^ (1ul << empty_spaces_down.at(i));

			std::vector<uShort> updated_empty_spaces_down = empty_spaces_down;
			//Remove the access for down electrons to be put here
			updated_empty_spaces_down.erase(updated_empty_spaces_down.begin(), updated_empty_spaces_down.begin() + i + 1);

			//Recursive
			combinationRecursiveAddingSingle(empty_spaces_up,updated_empty_spaces_down, added_num, n_up, n_down-1, sites, all_comb, placed);
		}

	}
}
void combinationAll(uShort n_up, uShort n_down, uShort sites, uShort nU, std::vector<sType>* all_states) {
	/****************************************************************
	* Computes the combinations of the states countaining <nUp> electrons up and <nDown> electrons down, while there is <nU> double occupation in a system of <sites> sites.
	* 
	* Parameter
	* ----------
	* n_up		: (uShort) number of up electrons to place
	* n_down	: (uShort) number of down electrons to place
	* sites		: (uShort) number of sites of the system
	* nU		: (uShort) number of double occupation
	* all_states: (vector<sType>*) Countains all the possible states
	*
	* Returns
	* -------
	* NONE
	*****************************************************************/
	//Number of states with specified nU, nUp, nDown values
	uLong combU_up_down = combSpecified(nU,sites,n_up,n_down);
	
	//Different ways to placed the double occupation
	uLong comb_U = comb(sites,nU);
	uLong comb_else = combU_up_down / comb_U;
	*all_states = std::vector<sType>(combU_up_down,0);

	if (n_up < nU || n_down < nU) {
		std::cout<<"Not enough electron up or down for the double occupation required"<<std::endl;
		exit(1);
	}

	if ((n_up + n_down - sites) > nU) return;

	//Independent electron number
	n_up -= nU;
	n_down -= nU;
	
	std::vector<sType> states_double_occupation;

	//Places the double occupation through the lattice
	if (nU > 0) {
		combinationDoubleOccupation(nU,sites,&states_double_occupation);
	}

	if (n_up > 0 || n_down > 0) {
		uLong placed = 0;
		#pragma omp parallel default(none) shared(states_double_occupation, all_states, sites, comb_else, n_up, n_down) private(placed)
		{
		#pragma omp for  
		for (uLong i = 0; i < states_double_occupation.size(); i++) {
			placed = 0;
			uLong current_num = states_double_occupation.at(i);

			//Finds where electron can still be placed
			std::vector<uShort> left_places;
			uLong one = 1; 
			for (uShort j = 0; j < sites; j++) {
				if((current_num & one) == 0) {
					left_places.push_back(j);
				}	
				one <<= 1;
			}
			//Add the leftover electron without double occupation
			combinationRecursiveAddingSingle(left_places, left_places, current_num,n_up, n_down, sites, all_states->data() + i * comb_else, &placed);

		}
		}

		if (nU == 0) {
			std::vector<uShort> chosen_sites(sites);
			std::iota(chosen_sites.begin(), chosen_sites.end(),0);
			uLong current_num = 0;
			placed = 0;
			combinationRecursiveAddingSingle(chosen_sites, chosen_sites,current_num,n_up,n_down,sites,all_states->data(),&placed);

		}
	}
	else { //If there is no loose electron
		for (uLong i = 0; i < states_double_occupation.size(); i++) {
			all_states->at(i) = (states_double_occupation.at(i));
		}
	}
}

bool acceptFunction(sType state, float acceptQuota){
	//std::srand(time(NULL));
	if (state == 0) std::cout<<"";//Line to shutup the warnings for now

	float rng_number = ((double)std::rand()) / RAND_MAX;

	bool accepted = (rng_number < acceptQuota);

	return accepted;
}

#pragma once
#include "binarySearchTree.h"


void isSuccess(bool status);
void cct(std::string text, int color);

int getNumberOfBlocks(int sites);

//Array search
double findMinOfArray(std::vector<double> array);

//Maths
unsigned long comb(uLong n, uLong r);
unsigned long comb_specified(uLong nU, uLong sites, uLong up, uLong down);
double fermi_dirac_function(float x, float beta,float mu);
double boltzmann_distribution_function(float dE, float beta);
double calculate_sd(std::vector<double> data);
int deg_fundamental_check(double* eigen_energies, uLong n,double eps = 10e-10);
template <class T, class R> void make_identity(T* matrix, R rows){
	/*****************************************
	* Creates the identity matrix out of the square matrix given
	*
	* Parameters
	* ----------
	* matrix: (T*) matrix to transform into identity.
	* rows	: (R) number of rows/columns of the square matrix
	*
	* Templates
	* ---------
	* T		: int, long, float, double, std::complex<double>
	* R		: int, long, unsigned 
	*
	* Returns
	* -------
	* NONE
	****************************************/
	for (R i = 0; i < rows; i++){
		for (R j = 0; j < rows; j++) {
			matrix[i * rows + j] = (i==j) ? 1 : 0;
		}
	}
}
template <class T, class U> void make_tri_diag(
    T* matrix, T* diag, T* h_diag, U size){
	/*****************************************
	* Creates a tridiagonal matrix with diag elements on the diagonal and h_diag ont the secondary diagonal
	* Parameters
	* ----------
	* matrix: (T*) matrix to transform into tridiag, must be given a matrix of zeros.
	* diag	: (T*) diagonal elements.
	* h_diag: (T*) secondary diagonal elements.
	* size	: (R) number of rows/columns of the square matrix
	*
	* Templates
	* ---------
	* T		: int, long, float, double, std::complex<double>
	* U		: int, long, unsigned 
	*
	* Returns
	* -------
	* NONE
	****************************************/
	for (U i = 0; i < size; i++){
		//diagonal elements
		matrix[i * (size + 1)] = diag[i];
		//Under diagonal
		if (i > 0){
			matrix[i*(size+1)-1] = h_diag[i-1];
		}
		//Over diagonal
		if (i<(size-1)){
			matrix[i*(size+1)+1] = h_diag[i];
		}
	}
}

//Vector manipulation
template <class R> void normalize(double* vec, R size) {
	/***********************************************************
	* Normalizes a vector

	* Parameters
	* ----------
	* vec	: (double*) vector to normalize
	* size	: (R) Number of elements of the vector
	*
	* Templates
	* ---------
	* R		: int, long, unsigned 
	* 
	* Returns
	* -------
	* NONE
	*************************************************************/
	//Normalisation of vec
	double norm = 1/cblas_dnrm2(size,vec,1);
	cblas_dscal(size,norm,vec,1);
}
template <class R> void initial_vector(
        R const SIZE, double* v, uInt seed=clock()) {
	/*********************************************************
	* Creates a random vector of size 'SIZE' and normalizes it
	*
	* Parameters
	* ----------
	* SIZE	: (R) vector dimension
	* v		: (double*) Receptacle for the random vector created
	* seed	: (unsigned int) Seed to generate random
	*
	* Templates:
	* ----------
	* T : float, double, std::complex<double>
	* R : unsigned int, unsigned long, ...
	*
	* Returns
	* -------
	* NONE
	*******************************************************/
	//Initialize seed
	srand(seed);
	for (R i = 0; i < SIZE; i++) {
		double a = (double)rand() / RAND_MAX;
		v[i] = a;
	}

	normalize(v,SIZE);
}

//Formating functions
std::string to_string_p(const double a_value, const int n=16);
template <class T> std::string to_string_pq(
    const T a_value, const uInt integers=18, const uInt decimals=12);

template <class T> T remove_zeros(T a);

//Print vec and matrix
template <class T, class R> std::string write_vector(const T* vec, R size){
	/*****************************************
	* Writes a double vector has a string 
	*
	* Parameters
	* ----------
	* vec	: (T*) pointer of the vector to stringify
	* size	: (R) size of the vector
	*
	* Templates:
	* ----------
	* T		: int, long, unsigned
	* R		: int, long, unsigned
	*
	* Returns
	* -------
	* write	: (std::string) double vector in string form
	****************************************/
	std::string write = "[";
	for (R i = 0; i < size; i++) {
		write += to_string_p(vec[i], 0) + "\t";
		if ((i+1) % 10 == 0){write += "\n";}
	}
	write += "]\n";
	return write;
}
template <class R> std::string write_matrix(
        const double* mat, R rows, R cols, 
        int integers = 2, int precision = 0){
	/*****************************************
	* Writes a double matrix has a string 
	*
	* Parameters
	* ----------
	* mat	    : (double*) pointer of the double matrix to stringify
	* rows	    : (R) number of rows of the matrix
	* cols	    : (R) number of cols of the matrix
    * integers  : (int) number of int spaces
    * precision : (int) number of float spaces
	*
	* Templates
	* ---------
	* R		: int, long, unsigned
	*
	* Returns
	* -------
	* write	: (std::string) double matrix in string form
	****************************************/
	std::string write = "";
	for (R i = 0; i < rows; i++) {
		write += "[";
		for (R j = 0; j < cols; j++) {
			write += to_string_pq(mat[i * cols + j],integers,precision)+"\t";
		}
		write += "]\n";
	}
	return write;
}

template <class T, class R> void print_vector(const T* vec, R size){
	/*****************************************
	* Prints a given vector 
	*
	* Parameters
	* ----------
	* mat	: (T*) pointer of the complex matrix to stringify
	* size	: (R) number of the vector
	*
	* Templates:
	* ----------
	* T		: double, std::complex<double>
	* R		: int, long, unsigned
	*
	* Returns
	* -------
	* NONE
	****************************************/
	std::cout << write_vector(vec, size);
}
template <class R> void print_matrix(const double* mat, R rows, R cols, 
                                     int integers = 2, int precision = 0){
	/*****************************************
	* Prints the given matrix
	*
	* Parameters
	* ----------
	* mat	: (double*) pointer of the complex matrix to stringify
	* rows	: (R) number of rows of the matrix
	* cols	: (R) number od cols of the matrix
    * integers  : (int) number of int spaces
    * precision : (int) number of float spaces
    *
	* Templates
	* ---------
	* T		: double, std::complex<double>
	* R		: int, long, unsigned
	*
	* Returns
	* -------
	* NONE
	****************************************/
	std::cout << write_matrix(mat, rows, cols, integers, precision);
}

std::string time_formating(
    std::chrono::time_point<std::chrono::system_clock> start, 
    std::chrono::time_point<std::chrono::system_clock> end);

unsigned int one_counter(sType num); 
//Bit operations
sType rotr(sType bitVec, unsigned char rotIndex, unsigned char size);



void combination_recursive(std::vector<uShort> empty_spaces, sType current_num, 
                           uShort left_to_place, sType* all_comb, sType placed, 
                           sType* it);

void combination_double_occupation(uShort N, uShort sites, 
                                   std::vector<sType>* all_double_occupation);
void combination_recursive_adding_single(
    std::vector<uShort> empty_spaces_up, std::vector<uShort> empty_spaces_down, 
    sType current_num, uShort n_up, uShort n_down, uShort sites, 
    sType* all_comb, sType* placed);
void combination_all(uShort n_up, uShort n_down, uShort sites, uShort nU, 
                     std::vector<sType>* all_states);

bool accept_function(sType state, float acceptQuota=0.5);

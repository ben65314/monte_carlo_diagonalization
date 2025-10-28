#ifndef __utilities_h__
#define __utilities_h__

#include <algorithm>
#include <bitset>
//#include <cblas.h>
#include <chrono>
#include <complex>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
//#include <lapack.h>
#include <limits.h>
#include <math.h>
#include <numeric>
#include <omp.h>
#include <random>
#include <sstream>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <tuple>
#include <type_traits>
#include <unistd.h>
#include <vector>
//
//For more than 16 sites
typedef uint64_t sType;		

extern "C" {

//BLAS
void dscal_(sType* n, double* alpha, double* x, int* incx);
//void dscal_(int, double, double*, int);
double dnrm2_(sType*, double*, int*);
void daxpy_(sType*, double*, double*, int*, double*, int*);
double ddot_(sType*, double*, int*, double*, int*);
void dswap_(sType*, double*, int*, double*, int*);

void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, 
          int*, double*, double*, int*);

//LAPACK
/* Subroutine */ int dsyev_(char *jobz, char *uplo, int *n, double *a,
	int *lda, double *w, double *work, int *lwork, 
	int *info);

/* Subroutine */ int dstev_(char *jobz, int *n, double *d__, 
	double *e, double *z__, int *ldz, double *work, 
	int *info);
}



///Types of the states;
//For 4 sites or less
//typedef uint8_t sType;

//For 8 sites or less
//typedef uint16_t sType;

//For 16 sites or less
//typedef uint32_t sType;


//Type of the vector 
typedef double vType;



typedef uint8_t uChar;
typedef uint16_t uShort;
typedef uint32_t uInt;
typedef uint64_t uLong;

//Declaration of external variables
extern int verbose;

#endif

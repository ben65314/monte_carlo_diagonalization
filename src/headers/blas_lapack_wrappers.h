#include "utilities.h"

extern "C" {

//BLAS definitions
//double
void dscal_(sType* n, double* alpha, double* x, int* incx);
double dnrm2_(sType*, double*, int*);
void daxpy_(sType*, double*, double*, int*, double*, int*);
double ddot_(sType*, double*, int*, double*, int*);
void dswap_(sType*, double*, int*, double*, int*);

void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*,
          int*, double*, double*, int*);
//complex
void zdscal_(sType* n, double* alpha, std::complex<double>* x, int* incx);
double dznrm2_(sType*, std::complex<double>*, int*);
void zaxpy_(sType*, std::complex<double>*, std::complex<double>*, int*, std::complex<double>*, int*);
double zdotc_(sType*, std::complex<double>*, int*, std::complex<double>*, int*);
void zswap_(sType*, std::complex<double>*, int*, std::complex<double>*, int*);
void zgemm_(char*, char*, int*, int*, int*, std::complex<double>*, std::complex<double>*, int*, std::complex<double>*,
          int*, std::complex<double>*, std::complex<double>*, int*);

//LAPACK
/* Subroutine */ int dsyev_(char *jobz, char *uplo, int *n, double *a,
	int *lda, double *w, double *work, int *lwork,
	int *info);

/* Subroutine */ int dstev_(char *jobz, int *n, double *d__,
	double *e, double *z__, int *ldz, double *work,
	int *info);

/* Subroutine */ int zheev_(char *jobz, char *uplo, int *n, std::complex<double>* a, int *lda, double* w, std::complex<double>* work, int* lwork, double* rwork, int *info);
}

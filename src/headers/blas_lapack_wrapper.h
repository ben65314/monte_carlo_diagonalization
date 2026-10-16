#include "utilities.h"
#ifndef __blas_lapack_wrapper_h__
#define __blas_lapack_wrapper_h__
extern "C" {

//BLAS definitions
//double

void dscal_(sType* n, double* alpha, double* x, int* incx);
void zdscal_(sType* n, double* alpha, std::complex<double>* x, int* incx);

double dnrm2_(sType*, double*, int*);
double dznrm2_(sType*, std::complex<double>*, int*);

void daxpy_(sType*, double*, double*, int*, double*, int*);
void zaxpy_(sType*, std::complex<double>*, std::complex<double>*, int*, std::complex<double>*, int*);

void dswap_(sType*, double*, int*, double*, int*);
void zswap_(sType*, std::complex<double>*, int*, std::complex<double>*, int*);

void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
void zgemm_(char*, char*, sType*, sType*, sType*, std::complex<double>*, std::complex<double>*, sType*, std::complex<double>*, sType*, std::complex<double>*, std::complex<double>*, sType*);

sType sdot_(sType*, sType*, int*, sType*, int*);
double ddot_(sType*, double*, int*, double*, int*);
std::complex<double> zdotc_(sType*, std::complex<double>*, int*, std::complex<double>*, int*);
void zdotcsub_(sType*, std::complex<double>*, int*, std::complex<double>*, int*, std::complex<double>* dotc);

//LAPACK
/* Subroutine */ int dstev_(char *jobz, int *n, double *d__,
	double *e, double *z__, int *ldz, double *work,
	int *info);
/* Subroutine */ int dsyev_(char *jobz, char *uplo, int *n, double *a,
	int *lda, double *w, double *work, int *lwork,
	int *info);
/* Subroutine */ int zheev_(char *jobz, char *uplo, int *n, std::complex<double>* a, int *lda, double* w, std::complex<double>* work, int* lwork, double* rwork, int *info);
}

//Wrappers
//Scaling of a vector
template <class D> void scal_blas(sType* n, double* scalar, D* vector) {
    int ONE = 1;
    if constexpr(std::is_same_v<D,double>) {
        dscal_(n, scalar, vector, &ONE);
    }
    else if constexpr(std::is_same_v<D,std::complex<double>>) {
        zdscal_(n, scalar, vector, &ONE);
    }
    else {static_assert(!sizeof(D), "scal_blas: unsupported type");}
}
//Norm of a vector
template <class D> double nrm2_blas(sType* n, D* vector) {
    int ONE = 1;
    double norm = 0;
    if constexpr(std::is_same_v<D,double>) {
        norm = dnrm2_(n, vector, &ONE);
    }
    else if constexpr(std::is_same_v<D,std::complex<double>>) {
        norm = dznrm2_(n, vector, &ONE);
    }
    else {static_assert(!sizeof(D), "nrm2_blas: unsupported type");}
    return norm;
}
//aX+Y
template <class D> void axpy_blas(sType* n, D* x, D* y, D a = 1) {
    int ONE = 1;
    if constexpr(std::is_same_v<D,double>) {
      daxpy_(n, &a, x, &ONE, y, &ONE);
    }
    else if constexpr(std::is_same_v<D,std::complex<double>>) {
        zaxpy_(n, &a, x, &ONE, y, &ONE);
    }
    else {static_assert(!sizeof(D), "axpy_blas: unsupported type");}
}
//Swap two vectors
template <class D> void swap_blas(sType* n, D* x, D* y) {
    int ONE = 1;
    if constexpr(std::is_same_v<D,double>) {
        dswap_(n, x, &ONE, y, &ONE);
    }
    else if constexpr(std::is_same_v<D,std::complex<double>>) {
        zswap_(n, x, &ONE, y, &ONE);
    }
    else {static_assert(!sizeof(D), "swap_blas: unsupported type");}
}
//Matrix product C = A*B
template <class D> void gemm_blas(char trans_a, char trans_b, int* m, int* n, int* k, D* A, D* B, D* C) {
    D a = 1; D b = 0;
    int* lda = (trans_a == 'N' || trans_a == 'n') ? m : k;
    int* ldb = (trans_a == 'N' || trans_a == 'n') ? k : n;
    int* ldc = m;
    if constexpr(std::is_same_v<D,double>) {
      dgemm_(&trans_a, &trans_b, m, n, k, &a, A, lda, B, ldb, &b, C, ldc);
    }
    else if constexpr(std::is_same_v<D,std::complex<double>>) {
        zgemm_(&trans_a, &trans_b, m, n, k, &a, A, lda, B, ldb, &b, C, ldc);
    }
    else {static_assert(!sizeof(D), "gemm_blas: unsupported type");}
}
//Dot product of two vectors
template <class D> D dot_blas(sType* n, D* x, D* y) {
    int ONE = 1;
    D dot = 0;
    if constexpr(std::is_same_v<D,double>) {
      dot = ddot_(n, x, &ONE, y, &ONE);
    }
    else if constexpr(std::is_same_v<D,std::complex<double>>) {
        zdotcsub_(n, x, &ONE, y, &ONE, &dot);
    }
    else {static_assert(!sizeof(D), "dot_blas: unsupported type");}
    return dot;
}
//Eigen solver
template <class D> void heev_lapack(char jobz, char uplo, int* n, D* A, double* w) {
    int lwork = *n * (*n+1);
    int info;
    D* work = new D[lwork];
    if constexpr(std::is_same_v<D,double>) {
      dsyev_(&jobz, &uplo, n, A, n, w, work, &lwork, &info);
    }
    else if constexpr(std::is_same_v<D,std::complex<double>>) {
        double* rwork = new double[3*(*n)-2];
        zheev_(&jobz, &uplo, n, A, n, w, work, &lwork, rwork, &info);
        delete[] rwork;
    }
    else {static_assert(!sizeof(D), "heev_lapack: unsupported type");}
    delete[] work;
}

#endif

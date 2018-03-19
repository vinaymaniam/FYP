
#ifndef ALGEBRA_H
#define ALGEBRA_H

#include <complex>
#include <stdexcept>
#include <vector>

#include "master.h"

#define REAL    1
#define COMPLEX 2

#ifndef _WIN32
#define F77_CALL(x) x ## _
#else
#define F77_CALL(x) x
#endif

#if MEX_COMPILE_FLAG==1
// include headers of Matlab's BLAS/LAPACK wrappers
// Matlab's wrapper use ptrdiff_t as integer type
typedef ptrdiff_t blasint_t;
#else
// declare extern Fortran wrappers for reference BLAS/LAPACK implementation
typedef int blasint_t;
#endif

#define XAXPY  CAT(PREFIX_REAL,F77_CALL(axpy))  // vector addition
#define XDOT   CAT(PREFIX_REAL,F77_CALL(dot))   // dot product
#define XGEMM  CAT(PREFIX_REAL,F77_CALL(gemm))  // matrix multiplication
#define ZGEMM  CAT(PREFIX_CPLX,F77_CALL(gemm))  // complex matrix multiplication
#define XSYMM  CAT(PREFIX_REAL,F77_CALL(symm))  // symmetric matrix multiplication
#define XSYRK  CAT(PREFIX_REAL,F77_CALL(syrk))  // multply matrix by its transpose
#define ZHERK  CAT(PREFIX_CPLX,F77_CALL(herk))  // multply complex matrix by its complex conjugate
#define XPOTRF CAT(PREFIX_REAL,F77_CALL(potrf)) // Choleski factorization of symmetric matrix
#define XPOTRI CAT(PREFIX_REAL,F77_CALL(potri)) // inverse of symmetric matrix
#define XPOTRS CAT(PREFIX_REAL,F77_CALL(potrs)) // solve linear symmetric system
#define XGEQRF CAT(PREFIX_REAL,F77_CALL(geqrf)) // QR decomposition  
#define XORMQR CAT(PREFIX_REAL,F77_CALL(ormqr)) // compute any of op(Q)*y or y*op(Q)
#define XTRTRS CAT(PREFIX_REAL,F77_CALL(trtrs)) // solve linear triangular system
#define XGESDD CAT(PREFIX_REAL,F77_CALL(gesdd)) // SVD decomposition
#define ZGESDD CAT(PREFIX_CPLX,F77_CALL(gesdd)) // complex SVD decomposition
#define XGEEV  CAT(PREFIX_REAL,F77_CALL(geev))  // computes eigenvalues and eigenvectors
#define ZGEEV  CAT(PREFIX_CPLX,F77_CALL(geev))  // computes complex eigenvalues and eigenvectors
#define XSYEVD CAT(PREFIX_REAL,F77_CALL(syevd)) // computes eigenvalues and eigenvectors for symmetric matrix
#define ZHEEVD CAT(PREFIX_CPLX,F77_CALL(heevd)) // computes compelx eigenvalues and eigenvectors for hermitian matrix

// BLAS/LAPACK function declaration (it allows to use any BLAS implementation
// without including external headers or change anything else in the code)
extern "C"{
void XAXPY(const blasint_t *n, const real_t *alpha, const real_t *x, 
        const blasint_t *incx, real_t *y, const blasint_t *incy);

real_t XDOT(const blasint_t *n, const real_t *dx, const blasint_t *incx, 
        const real_t *dy, const blasint_t *incy);
    
void XGEMM(char *transa, char *transb, const blasint_t *n, const blasint_t *m, 
        const blasint_t *k, real_t *alpha, const real_t *A, const blasint_t *lda, const real_t *B, 
        const blasint_t *ldb, real_t *beta, real_t *C, const blasint_t *ldc);

void ZGEMM(char *transa, char *transb, const blasint_t *n, const blasint_t *m, 
        const blasint_t *k, real_t *alpha, const real_t *A, const blasint_t *lda, const real_t *B, 
        const blasint_t *ldb, real_t *beta, real_t *C, const blasint_t *ldc);

void XSYMM(char *side, char *uplo, const blasint_t *m, const blasint_t *n, 
        real_t *alpha, const real_t *A, const blasint_t *lda, const real_t *B, 
        const blasint_t *ldb, real_t *beta, real_t *C, const blasint_t *ldc);

void XSYRK(char *uplo, char *trans, const blasint_t *n, const blasint_t *k, 
        real_t *alpha, const real_t *A, const blasint_t *lda, real_t *beta, 
        real_t *C, const blasint_t *ldc);

void ZHERK(char *uplo, char *trans, const blasint_t *n, const blasint_t *k, 
        real_t *alpha, const real_t *A, const blasint_t *lda, real_t *beta, 
        real_t *C, const blasint_t *ldc);

void XPOTRF(char *uplo, const blasint_t *n, real_t *A, const blasint_t *lda, 
        blasint_t *info);

void XPOTRS(char *uplo, const blasint_t *n, const blasint_t *nrhs, 
        const real_t *A, const blasint_t *lda, real_t *B, const blasint_t *ldb,
        blasint_t *info);

void XPOTRI(char *uplo, const blasint_t *n, real_t *A, const blasint_t *lda, 
        blasint_t *info);

void XGEQRF(const blasint_t *m, const blasint_t *n, real_t *A, const blasint_t *lda, 
        real_t *tau, real_t *work, const blasint_t *lwork, blasint_t *info);

void XORMQR(char *side, char *trans, const blasint_t *m, const blasint_t *n, 
        const blasint_t *k, const real_t *A, const blasint_t *lda, const real_t *tau, real_t *C, 
        const blasint_t *ldc, real_t *work, const blasint_t *lwork, blasint_t *info);

void XTRTRS(char *uplo, char *trans, char *diag, const blasint_t *n, 
        const blasint_t *nrhs, const real_t *A, const blasint_t *lda, real_t *B, const blasint_t *ldb, 
        blasint_t *info);

void XGESDD(const char *jobz, const blasint_t *m, const blasint_t *n, real_t *A, 
        const blasint_t *lda, real_t *S, real_t *U, const blasint_t *ldu, real_t *Vt, 
        const blasint_t *ldvt, real_t *work, const blasint_t *lwork, const blasint_t *iwork, 
        blasint_t *info);

void ZGESDD(const char *jobz, const blasint_t *m, const blasint_t *n, real_t *A, 
        const blasint_t *lda, real_t *S, real_t *U, const blasint_t *ldu, real_t *Vt, 
        const blasint_t *ldvt, real_t *work, const blasint_t *lwork, real_t *rwork, 
        const blasint_t *iwork, blasint_t *info);

void XGEEV(const char *jobvl, const char *jobvr, const blasint_t *n, real_t *A, 
        const blasint_t *lda, real_t *wr, real_t *wi, real_t *VL, const blasint_t *ldvl, 
        real_t *VR, const blasint_t *ldvr, real_t *work, const blasint_t *lwork, 
        blasint_t *info);

void ZGEEV(const char *jobvl, const char *jobvr, const blasint_t *n, const real_t *A, 
        const blasint_t *lda, real_t *w, real_t *VL, const blasint_t *ldvl, 
        real_t *VR, const blasint_t *ldvr, real_t *work, const blasint_t *lwork, 
        real_t *rwork, blasint_t *info);

void XSYEVD(const char *jobz, const char *uplo, const blasint_t *n, real_t *A, 
        const blasint_t *lda, real_t *W, real_t *work, const blasint_t *lwork, 
        blasint_t *iwork, const blasint_t *liwork, blasint_t *info);

void ZHEEVD(const char *jobz, const char *uplo, const blasint_t *n, real_t *A, 
        const blasint_t *lda, real_t *W, real_t *work, const blasint_t *lwork, 
        real_t *rwork, const blasint_t *lrwork, blasint_t *iwork, 
        const blasint_t *liwork, blasint_t *info);
}

/* NOTES:
   BLAS/LAPACK stores matrix in vectorized column-major form
   an array of complex numbers in LAPACK is an array of double size with the 
   real and imaginary parts of each complex number interleaved */

/* transform std::complex vector x = {(r1,i1),(r2,i2),...,(rn,in)} to its 
   fortran representation y = {r1,i1,r2,i2,...,rn,in}  */
std::vector<real_t> c2fort(const std::vector<complex_t> &x) {
    std::vector<real_t> y;
    y.reserve(2*x.size());
    for(size_t i=0; i<x.size(); i++) {
        y.push_back( x[i].real() );
        y.push_back( x[i].imag() );
    }
    return y;
}

/* transform a vector x = {r1,i1,r2,i2,...,rn,in} from its complex fortran 
   representation to std::complex y = {(r1,i1),(r2,i2),...,(rn,in)} */
std::vector<complex_t> fort2c(const std::vector<real_t> &x) {
    #if SAFE_CHECK
    if(x.size()%2!=0) {
        throw std::invalid_argument("fort2c : invalid dimensions");
    }
    #endif
    std::vector<complex_t> y;
    y.reserve(x.size()/2);
    for(size_t i=0; i<x.size(); i+=2) {
        y.push_back( {x[i], x[i+1]} );
    }
    return y;
}

/* calculates vector operation alpha*x+y; y is overwritten with the result */
void axpy(size_t n, const real_t alpha, const std::vector<real_t> &x, std::vector<real_t> &y) {
    #if SAFE_CHECK
    if(y.size()<x.size()) {
        throw std::invalid_argument("axpy : invalid output dimension");
    }
    #endif
    const blasint_t nn = (blasint_t)n;
    const blasint_t inc = 1;
    XAXPY(&nn, &alpha, &x[0], &inc, &y[0], &inc);
}

/* calculates dot product <x*y> */
real_t dot(size_t n, const std::vector<real_t> &x, const std::vector<real_t> &y) {
    #if SAFE_CHECK
    if(y.size()!=x.size()) {
        throw std::invalid_argument("dot : invalid input dimension");
    }
    #endif
    const blasint_t nn = (blasint_t)n;
    const blasint_t inc = 1;
    return XDOT(&nn, &x[0], &inc, &y[0], &inc);
}

/* calculate matrix product C = A*B where op(A) is m-by-k, op(B) is k-by-n, 
   and C is m-by-n, being op(X)=X if transx='N' and op(X)=X' if transx='T';
   matrices A, B, C are stored in memory in column-major vectorized form,
   assume A is M-by-N matrix vectorized in memory as an array of size M*N:
   A(r1:M-r2,c1:N-c2) => m=M-r2-r1+1  lda=M   k=N-c2-c1+1   osa=c1*M+r1 */
void gemm(int mode, char transa, char transb, size_t m, size_t n, size_t k, const std::vector<real_t> &A, 
        size_t lda, const std::vector<real_t> &B, size_t ldb, std::vector<real_t> &C, size_t ldc, 
        size_t osa = 0, size_t osb = 0, size_t osc = 0) {
    #if SAFE_CHECK
    if(C.size()<m*n*mode) {
        throw std::invalid_argument("gemm : invalid output dimension");
    }
    #endif
    std::vector<real_t> zero = {0.0, 0.0};
    std::vector<real_t> one = {1.0, 0.0};
    const blasint_t mm = (blasint_t)m;
    const blasint_t nn = (blasint_t)n;
    const blasint_t kk = (blasint_t)k;
    const blasint_t llda = (blasint_t)lda;
    const blasint_t lldb = (blasint_t)ldb;
    const blasint_t lldc = (blasint_t)ldc;
    if(mode==COMPLEX) {
        ZGEMM(&transa, &transb, &mm, &nn, &kk, &one[0], 
                &A[osa], &llda, &B[osb], &lldb, &zero[0], &C[osc], &lldc);
    } else {
        XGEMM(&transa, &transb, &mm, &nn, &kk, &one[0], 
                &A[osa], &llda, &B[osb], &lldb, &zero[0], &C[osc], &lldc);
    }
}

/* calculate matrix product C=A*B (when side='L') or C=B*A (when size='R')
   being C a m-by-n matrix, A a symmetric m-by-m (when side='L') or n-by-n
   (when side='R') matrix with values filled in the upper triangular part
   (when uplo='U') or lower triangular (when uplo='L') */
void symm(char side, char uplo, size_t m, size_t n, const std::vector<real_t> &A, size_t lda, 
        const std::vector<real_t> &B, size_t ldb, std::vector<real_t> &C, size_t ldc) {
    #if SAFE_CHECK
    if(C.size()<m*n) {
        throw std::invalid_argument("symm : invalid output dimension");
    }
    #endif
    real_t zero = 0.0;
    real_t one = 1.0;
    const blasint_t mm = (blasint_t)m;
    const blasint_t nn = (blasint_t)n;
    const blasint_t llda = (blasint_t)lda;
    const blasint_t lldb = (blasint_t)ldb;
    const blasint_t lldc = (blasint_t)ldc;
    XSYMM(&side, &uplo, &mm, &nn, &one, &A[0], &llda, &B[0], &lldb, &zero, &C[0], &lldc);
}

/* computes C=A*A' (when transa='N') or C=A'*A (when transa='T'), being C 
   a n-by-n matrix,and A a n-by-k (when transa='N') or k-by-n (when transa='T')
   matrix; if A is complex, A' will be its complex conjugate; on exit the 
   upper (when uplo='U') or lower (when uplo='L') triangular part of C will 
   be filled with the result */
void syrk(int mode, char uplo, char transa, size_t n, size_t k, const std::vector<real_t> &A, 
        size_t lda, std::vector<real_t> &C, size_t ldc) {
    #if SAFE_CHECK
    if(C.size()<n*n*mode) {
        throw std::invalid_argument("dsyrk : invalid output dimension");
    }
    #endif
    std::vector<real_t> zero = {0.0, 0.0};
    std::vector<real_t> one = {1.0, 0.0};
    const blasint_t nn = (blasint_t)n;
    const blasint_t kk = (blasint_t)k;
    const blasint_t llda = (blasint_t)lda;
    const blasint_t lldc = (blasint_t)ldc;
    if(mode==COMPLEX) {
        ZHERK(&uplo, &transa, &nn, &kk, &one[0], &A[0], &llda, &zero[0], &C[0], &lldc);
    } else {
        XSYRK(&uplo, &transa, &nn, &kk, &one[0], &A[0], &llda, &zero[0], &C[0], &lldc);
    }
}

/* computes the Cholesky factorization of a real symmetric positive definite 
   matrix A of size n-by-n with values filled in the upper triangular part
   (when uplo='U') or lower triangular (when uplo='L'); on exit the
   triangular part defined by uplo is overwritten with result U or L which
   satisfy either A = U'*T*U or A = L*L'*T, being U and L the upper or
   lower triangular factor matrices of A, respectively */
int potrf(char uplo, size_t n, std::vector<real_t> &A, size_t lda) {
    #if SAFE_CHECK
    if(A.size()<n*n) {
        throw std::invalid_argument("potrf : invalid input dimension");
    }
    #endif
    const blasint_t nn = (blasint_t)n;
    const blasint_t llda = (blasint_t)lda;
    blasint_t info = 0;
    XPOTRF(&uplo, &nn, &A[0], &llda, &info);
    #if SAFE_CHECK
    if(info<0) {
        throw std::runtime_error("potrf : illegal value");
    }
    if(info>0) {
        throw std::runtime_error("potrf : matrix is non positive definite");
    }
    #endif
    return (int)info;
}

/* compute the inverse of a real symmetric positive definite matrix n-by-n,
   whose Cholesky factorization U or L, as computed by potrf, is stored in 
   the upper triangular part (when uplo='U') or lower triangular (when 
   uplo='L') of A; on exit the upper or lower triangular part of A (i.e., 
   the input factor) is overwritten with the result (i.e., the symmetric 
   inverse matrix) */
int potri(char uplo, size_t n, std::vector<real_t> &A, size_t lda) {
    #if SAFE_CHECK
    if(A.size()<n*n) {
        throw std::invalid_argument("potri : invalid input dimension");
    }
    #endif
    const blasint_t nn = (blasint_t)n;
    const blasint_t llda = (blasint_t)lda;
    blasint_t info = 0;
    XPOTRI(&uplo, &nn, &A[0], &llda, &info);
    #if SAFE_CHECK
    if(info<0) {
        throw std::runtime_error("potri : illegal value");
    }
    if(info>0) {
        throw std::runtime_error("potri : zero element in factorization");
    }
    #endif
    return (int)info;
}

/* solve a system of linear equations A*X = B with a symmetric positive 
   definite matrix A n-by-n using the Cholesky factorization A = U'*T*U or 
   A = L*L'*T computed by potrf; the upper (when uplo='U') or lower (when
   uplo='L') triangular part of A contains the triangular factor U or L 
   from the Cholesky decomposition, nrhs is the number of columns in B (i.e.
   the number of right hand sides in the linear equation), B is the right 
   hand side of the equation; on exit B is overwritten with the result X */
int potrs(char uplo, size_t n, size_t nrhs, const std::vector<real_t> &A, 
        size_t lda, std::vector<real_t> &B, size_t ldb) {
    #if SAFE_CHECK
    if(B.size()<n*n) {
        throw std::invalid_argument("potrs : invalid input dimension");
    }
    #endif
    const blasint_t nn = (blasint_t)n;
    const blasint_t nnrhs = (blasint_t)nrhs;
    const blasint_t llda = (blasint_t)lda;
    const blasint_t lldb = (blasint_t)ldb;
    blasint_t info = 0;
    XPOTRS(&uplo, &nn, &nnrhs, &A[0], &llda, &B[0], &lldb, &info);
    #if SAFE_CHECK
    if(info<0) {
        throw std::runtime_error("potrs : illegal value");
    }
    #endif
    return (int)info;
}

/* compute QR decomposition of matrix A; on exit, the elements on and above
   the diagonal of the array A contain the min(m,n)-by-n upper trapezoidal 
   matrix R (R is upper triangular if m >= n); the elements below the 
   diagonal, with the array TAU, represent the orthogonal matrix Q as a
   product of min(m,n) elementary reflectors */
int geqrf(size_t m, size_t n, std::vector<real_t> &A, size_t lda, std::vector<real_t> &tau) {
    #if SAFE_CHECK
    if(tau.size()<n) {
        throw std::invalid_argument("geqrf : invalid output dimension");
    }
    #endif
    const blasint_t mm = (blasint_t)m;
    const blasint_t nn = (blasint_t)n;
    const blasint_t llda = (blasint_t)lda;
    blasint_t info = 0;
    // query optimal workspace size
    real_t size;
    const blasint_t query = -1;
    XGEQRF(&mm, &nn, &A[0], &llda, &tau[0], &size, &query, &info);
    // calculate QR
    const blasint_t lwork = (blasint_t)size;
    std::vector<real_t> work(lwork);
    XGEQRF(&mm, &nn, &A[0], &llda, &tau[0], &work[0], &lwork, &info);
    #if SAFE_CHECK
    if(info<0) {
        throw std::runtime_error("geqrf : illegal value");
    }
    #endif
    return (int)info;
}

/* compute any of Q'*C, Q*C, C*Q, or C*Q', this can be used to solve linear 
   least-square problem being Q orthogonal, Q is represented as product of 
   reflectors (as returned by XGEQRF), the ith column of A must contain the 
   ith elementary reflector, on exit C is overwritten by the solution, 
   A is m-by-k with lda, and C is m-by-n with ldc; parameter side constrols
   if Q is applied to the left or to the right, and trans controls if Q
   is transposed ('T') or not ('N') */
int ormqr(char side, char trans, size_t m, size_t n, size_t k, const std::vector<real_t> &A, size_t lda,
        const std::vector<real_t> &tau, std::vector<real_t> &C, size_t ldc) {
    #if SAFE_CHECK
    if(C.size()<n) {
        throw std::invalid_argument("ormqr : invalid output dimension");
    }
    #endif
    const blasint_t mm = (blasint_t)m;
    const blasint_t nn = (blasint_t)n;
    const blasint_t kk = (blasint_t)k;
    const blasint_t llda = (blasint_t)lda;
    const blasint_t lldc = (blasint_t)ldc;
    blasint_t info = 0;
    // query optimal workspace size
    real_t size;
    const blasint_t query = -1;
    XORMQR(&side, &trans, &mm, &nn, &kk, &A[0], &llda, 
            &tau[0], &C[0], &lldc, &size, &query, &info);
    // calculate product
    const blasint_t lwork = (blasint_t)size;
    std::vector<real_t> work(lwork);
    XORMQR(&side, &trans, &mm, &nn, &kk, &A[0], &llda, 
            &tau[0], &C[0], &lldc, &work[0], &lwork, &info);
    #if SAFE_CHECK
    if(info<0) {
        throw std::runtime_error("ormqr : illegal value");
    }
    #endif
    return (int)info;
}

/* solve a triangular system of the form Ax=B, A is an upper (uplo='U') or
   lower (uplo='L') triangular matrix of size n-by-n, parameters trans 
   controls whether A is transposed ('T') or not ('N'), and diag controls 
   if A is non-unit ('N') or unit ('U') diagonal; B is n-by-nrhs, and on 
   exit it contains the solution */
int trtrs(char uplo, char trans, char diag, size_t n, size_t nrhs, 
        const std::vector<real_t> &A, size_t lda, std::vector<real_t> &B, size_t ldb) {
    #if SAFE_CHECK
    if(B.size()<n) {
        throw std::invalid_argument("trtrs : invalid output dimension");
    }
    #endif
    const blasint_t nn = (blasint_t)n;
    const blasint_t nnrhs = (blasint_t)nrhs;
    const blasint_t llda = (blasint_t)lda;
    const blasint_t lldb = (blasint_t)ldb;
    blasint_t info = 0;
    XTRTRS(&uplo, &trans, &diag, &nn, &nnrhs, &A[0], &llda, &B[0], &lldb, &info);
    #if SAFE_CHECK
    if(info<0) {
        throw std::runtime_error("trtrs : illegal value");
    }
    if(info>0) {
        throw std::runtime_error("trtrs : singular matrix");
    }
    #endif
    return (int)info;
}

/* singular value decomposition A = U*S*V', where A is m-by-n, U is m-by-m,
   S is a real vector 1-by-min(m,n) and V' is n-by-n, S is always real whereas
   U and V' have the same numeric class as A (real or complex), A is not passed
   by reference and it is not a constant because it is overwritten by xGESDD;
   if jobz='A', then all m columns of U and all n rows of Vt are computed, 
   if jobz='S', then the first min(m,n) columns of U and rows of Vt are computed*/
int gesdd(int mode, char jobz, size_t m, size_t n, std::vector<real_t> A, size_t lda,
        std::vector<real_t> &S, std::vector<real_t> &U, std::vector<real_t> &Vt) {
    #if SAFE_CHECK
    if(S.size()<MIN(m,n) || U.size()<m*m*mode || Vt.size()<n*n*mode) {
        throw std::invalid_argument("gesdd : invalid output dimension");
    }
    #endif
    if(mode!=REAL && mode!=COMPLEX) {
        mode = REAL;
    }
    const size_t k = MIN(m,n);
    const blasint_t mm = (blasint_t)m;
    const blasint_t nn = (blasint_t)n;
    const blasint_t llda = (blasint_t)lda;
    std::vector<real_t> rwork(5*k*(k+1)); // from LAPACK definition (maybe 7*k)
    std::vector<blasint_t> iwork(8*k); // from LAPACK definition
    blasint_t info = 0;
    // query optimal workspace size
    real_t size;
    const blasint_t query = -1;
    if(mode==COMPLEX) {
        ZGESDD(&jobz, &mm, &nn, &A[0], &llda, &S[0], &U[0], &mm, &Vt[0], &nn, 
                &size, &query, &rwork[0], &iwork[0], &info);
    } else {
        XGESDD(&jobz, &mm, &nn, &A[0], &llda, &S[0], &U[0], &mm, &Vt[0], &nn, 
                &size, &query, &iwork[0], &info);
    }
    // calculate SVD
    const blasint_t lwork = (blasint_t)size;
    std::vector<real_t> work(lwork*mode);
    if(mode==COMPLEX) {
        ZGESDD(&jobz, &mm, &nn, &A[0], &llda, &S[0], &U[0], &mm, &Vt[0], &nn, 
                &work[0], &lwork, &rwork[0], &iwork[0], &info);
    } else {
        XGESDD(&jobz, &mm, &nn, &A[0], &llda, &S[0], &U[0], &mm, &Vt[0], &nn, 
                &work[0], &lwork, &iwork[0], &info);
    }
    #if SAFE_CHECK
    if(info<0) {
        throw std::runtime_error("gesdd : illegal value");
    }
    if(info>0) {
        throw std::runtime_error("gesdd : algorithm could not converge");
    }
    #endif
    return (int)info;
}

/* compute pseudoinverse X n-by-m of A m-by-n such that X has the same 
   dimensions of A' and A*X*A = A, X*A*X = X and A*X and X*A are Hermitian */
int pinv(int mode, size_t m, size_t n, const std::vector<real_t> &A, size_t lda, std::vector<real_t> &X) {
    #if SAFE_CHECK
    if(X.size()<m*n*mode) {
        throw std::invalid_argument("pinv : invalid output dimension");
    }
    #endif
    const size_t k = MIN(m,n);
    std::vector<real_t> S(k);
    std::vector<real_t> U(m*m*mode);
    std::vector<real_t> Vt(n*n*mode);
    // calculate SVD of A
    int info = gesdd(mode, 'S', m, n, A, lda, S, U, Vt);
    // equivalent to Matlab eps(max(S)), XGESDD returns singular values
    // sorted in descending order, thus max(S)=S[0]
    real_t factor = EPSILON * (real_t)(1<<(int)std::floor(std::log2( S[0] )));
    // ignore all eigenvalues smaller than tol
    real_t tol = MAX(m,n) * factor;
    size_t r = 1;
    // only first k rows (right singular vectors) of Vt are filled in 'S' mode
    // out of these k, we select the rows corresponding to singular values>=tol 
    // we denote the number of rows with singular values>=tol as r
    // Vt is vectorized column-major, thus Vt is n-by-n with ldv=r (being r<=k)
    for(size_t i=1; i<k && S[i]>=tol; i++) {
        r++;
    }
    #if SAFE_CHECK
    if(r==0 || r>k) {
        throw std::invalid_argument("pinv : invalid singular values number : tol=" + 
                std::to_string(tol) + " / max(S)=" + std::to_string(S[0]));
    }
    #endif
    // NOTE: r must be at least 1
    for(size_t i=0; i<r*mode; i+=mode) {
        for(size_t j=0; j<n; j++) {
            // divide each element j in row i by ith singular value 
            Vt[j*n*mode+i] = Vt[j*n*mode+i] / S[i/mode];
            if(mode==COMPLEX) {
                Vt[j*n*mode+i+1] = Vt[j*n*mode+i+1] / S[i/mode];
            }
        }
    }
    // calculate X = Vt(1:r,:)'*U(:,1:r)' where 
    //  Vt(1:r,:)' is n-by-r with lda=n
    //  U(:,1:r)' is r-by-m with ldb=m
    //  X is the pseudo-inverse n-by-m with ldc=n (vectorized column-major)
    //  conjugate transposition if data is complex-valued
    char trans = (mode==COMPLEX)? 'C' : 'T';
    gemm(mode, trans, trans, n, m, r, Vt, n, U, m, X, n);
    return info;
}

/* compute (complex) eigenvalues W n-by-1 of non-symmetric real or complex  
   matrix A n-by-n; if jobv is 'N' left/right eigenvectors are not computed,
   whereas if jobv is 'V' left/right eigenvectors are computed; A is 
   overwritten on exit with the result */
int geev(int mode, char jobv, size_t n, std::vector<real_t> &A, size_t lda, std::vector<real_t> &W) {
    #if SAFE_CHECK
    if(W.size()<2*n) {
        throw std::invalid_argument("geev : invalid output dimension");
    }
    #endif
    const blasint_t nn = (blasint_t)n;
    const blasint_t llda = (blasint_t)lda;
    std::vector<real_t> WR; // contains real part of eigenvalues 
    std::vector<real_t> WI; // contains imaginary part of eigenvalues
    std::vector<real_t> VLR(1); // eigenvectors (non-referenced if jobv='N')
    const blasint_t ldv = 1;
    std::vector<real_t> rwork;
    blasint_t info = 0;
    // query optimal workspace size
    std::vector<real_t> size(2);
    const blasint_t query = -1;
    if(mode==COMPLEX) {
        rwork = std::vector<real_t>(2*n); // from LAPACK definition
        ZGEEV(&jobv, &jobv, &nn, &A[0], &llda, &W[0], &VLR[0], &ldv, &VLR[0], &ldv,
                &size[0], &query, &rwork[0], &info);
    } else {
        WR = std::vector<real_t>(n);
        WI = std::vector<real_t>(n);
        XGEEV(&jobv, &jobv, &nn, &A[0], &llda, &WR[0], &WI[0], &VLR[0], &ldv, 
                &VLR[0], &ldv, &size[0], &query, &info);
    }
    // caclulate eigenvalues
    const blasint_t lwork = (blasint_t)size[0];
    std::vector<real_t> work(lwork*mode);
    if(mode==COMPLEX) {
        ZGEEV(&jobv, &jobv, &nn, &A[0], &llda, &W[0], &VLR[0], &ldv, &VLR[0], &ldv,
                &work[0], &lwork, &rwork[0], &info);
    } else {
        XGEEV(&jobv, &jobv, &nn, &A[0], &llda, &WR[0], &WI[0], &VLR[0], &ldv, 
                &VLR[0], &ldv, &work[0], &lwork, &info);
        // merge real and imaginary parts in output vector W
        for(size_t i=0; i<n; i++) {
            W[2*i]   = WR[i];
            W[2*i+1] = WI[i];
        }
    }
    #if SAFE_CHECK
    if(info<0) {
        throw std::runtime_error("geev : illegal value");
    }
    if(info>0) {
        throw std::runtime_error("geev : algorithm could not converge");
    }
    #endif
    return (int)info;
}

/* compute real-valued eigenvalues W n-by-1 of real symmetric or complex 
   Hermitian matrix A n-by-n with values filled in the upper triangular part
   (when uplo='U') or lower triangular (when uplo='L'); on exit, if jobv='V'
   A is overwritten with the complex orthogonal eigenvectors (stored 
   column-wise); if jobv='N' the content of A is destroyed; W on exit will 
   contain the eigenvalues in ascending order */
int syevd(int mode, char jobv, char uplo, size_t n, std::vector<real_t> &A, 
        size_t lda, std::vector<real_t> &W) {
    #if SAFE_CHECK
    if(W.size()<n) {
        throw std::invalid_argument("syevd : invalid output dimension");
    }
    #endif
    const blasint_t nn = (blasint_t)n;
    const blasint_t llda = (blasint_t)lda;
    blasint_t iwkopt;
    real_t rwkopt = -1;
    std::vector<real_t> wkopt(mode);
    blasint_t info = 0;
    // query optimal workspace size
    blasint_t lwork = -1;
    blasint_t lrwork = -1;
    blasint_t liwork = -1;
    if(mode==COMPLEX) {
        // the routine only calculates the optimal sizes of the WORK, RWORK, 
        // and IWORK arrays, which are returned in the first entries of the 
        // WORK, RWORK and IWORK arrays
        ZHEEVD(&jobv, &uplo, &nn, &A[0], &llda, &W[0], &wkopt[0], &lwork, &rwkopt,
                &lrwork, &iwkopt, &liwork, &info);
    } else {
        // the routine only calculates the optimal sizes of the WORK and 
        // IWORK arrays, which are returned in the first entries of the 
        // WORK and IWORK arrays
        XSYEVD(&jobv, &uplo, &nn, &A[0], &llda, &W[0], &wkopt[0], &lwork, 
                &iwkopt, &liwork, &info);
    }
    // solve eigenproblem
    lwork = (blasint_t)wkopt[0];
    std::vector<real_t> work(lwork*mode);
    lrwork = (blasint_t)rwkopt;
    std::vector<real_t> rwork(lrwork);
    liwork = iwkopt;
    std::vector<blasint_t> iwork(liwork);
    if(mode==COMPLEX) {
        ZHEEVD(&jobv, &uplo, &nn, &A[0], &llda, &W[0], &work[0], &lwork, 
                &rwork[0], &lrwork, &iwork[0], &liwork, &info);
    } else {
        XSYEVD(&jobv, &uplo, &nn, &A[0], &llda, &W[0], &work[0], &lwork, 
                &iwork[0], &liwork, &info);
    }
    #if SAFE_CHECK
    if(info<0) {
        throw std::runtime_error("syevd : illegal value");
    }
    if(info>0) {
        throw std::runtime_error("syevd : algorithm could not converge");
    }
    #endif
    return (int)info;
}

#endif

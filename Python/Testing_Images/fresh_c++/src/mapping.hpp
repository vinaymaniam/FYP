
#ifndef MAPPING_H
#define MAPPING_H

#include <stdexcept>
#include <vector>

#include "master.h"
#include "algebra.hpp"

// modifiable constants
#define LIN_MAP_MODE 3  // integer {1->pinv, 2->Cholesky inv, 3->symm system}

/* apply linear mapping; the goal is to calculate the linear mapping matrix 
   M, i.e. the model for the local reconstruction learned from the local 
   training sets Yh and Yl (containig low- and high-res mutually similar 
   example patches); the patches can have arbitrary dimension, and are 
   assumed to be given in vectorized form in a vector s-by-1; the model M 
   is used to predict the high-res version of the reference patch as
                   M = (Yh*Yl') * (Yl*Yl'+lambda*I)^-1     (1)
                     =      Yhl * Yll^-1                   (2)
   where Yh, Yl are s-by-n matrices (each column contains a vectorized 
   patch of length s), I is the identity matrix s-by-s and lambda>=0 is a 
   regularization parameter; both Yhl and Yll are s-by-s matrices, thus M 
   is also s-by-s. NOTE: the inversion of the matrix Yll is critical */
inline std::vector<real_t> linear_mapping(const std::vector<real_t> &xl, 
        const std::vector<real_t> &Yl, const std::vector<real_t> &Yh, 
        size_t s, size_t n, real_t lambda) {
    
    #if SAFE_CHECK
    if(Yl.size()!=s*n || Yh.size()!=s*n || xl.size()!=s) {
        throw std::invalid_argument("linear_mapping : invalid data size");
    }
    if(Yl.size()!=Yh.size()) {
        throw std::invalid_argument("linear_mapping : group dimensions do not agree");
    }
    #endif

    // pre-allocation
    std::vector<real_t> xh(s);
    std::vector<real_t> Yhl(SQUARE(s));
    std::vector<real_t> Yll(SQUARE(s));
    #if LIN_MAP_MODE!=3
    std::vector<real_t> M(SQUARE(s));
    #endif

    // direct (VERY inefficient) derivation of M using pinv
    #if LIN_MAP_MODE==1
    // Yl*Yl'
    gemm(REAL, 'N', 'T', s, s, n, Yl, s, Yl, s, Yhl, s);
    // Yl*Yl'+lambda*I
    for(size_t d=0; d<s; d++) {
        Yhl[d*(s+1)] += lambda;
    }
    // Yll <- pinv(Yl*Yl'+lambda*I) [pseudo-inverse]
    pinv(REAL, s, s, Yhl, s, Yll);
    // Yhl <- Yh*Yl'
    gemm(REAL, 'N', 'T', s, s, n, Yh, s, Yl, s, Yhl, s);
    // M <- Yhl*Yll
    gemm(REAL, 'N', 'N', s, s, s, Yhl, s, Yll, s, M, s);
    // xh <- M*xl
    gemm(REAL, 'N', 'N', s, 1, s, M, s, xl, s, xh, s);

    // let us note that Yll=Yl*Yl' is (by definition) symmetric and 
    // positive definite, thus there is a more effcient way to compute 
    // Yll^-1 using Cholesky decomposition on Yll; moreover we can exploit 
    // several LAPACK routine specifically designed for symmetric 
    // input/output which only address the upper (or lower) triangular 
    // part of the symmetric matrix effectively reducing the required 
    // computations by a factor of 2
    #elif LIN_MAP_MODE==2
    // Yl*Yl'
    syrk(REAL, 'U', 'N', s, n, Yl, s,  Yll, s);
    // Yl*Yl'+lambda*I
    for(size_t d=0; d<s; d++) {
        Yll[d*(s+1)] += lambda;
    }
    // Cholesky factorization of Yl*Yl'+lambda*I
    potrf('U', s, Yll, s);
    // Yll <- (Yl*Yl'+lambda*I)^-1 [inverse Cholesky]
    potri('U', s, Yll, s);
    // Yhl <- Yh*Yl'
    gemm(REAL, 'N', 'T', s, s, n, Yh, s, Yl, s, Yhl, s);
    // M <- Yhl*Yll 
    symm('R', 'U', s, s, Yll, s, Yhl, s, M, s);
    // xh <- M*xl
    gemm(REAL, 'N', 'N', s, 1, s, M, s, xl, s, xh, s);

    // in general, the inversion of Yll might be numerically unstable 
    // (especially in single/float precision); thus it is good practice to 
    // always avoid direct matrix inversion; instead let us solve the 
    // underlying linear system, where M is the unknown: thus the system 
    // (2) has a non-canonical form XA=B; we can rewrite it in a canonical 
    // form AX=B by transpoing both sides (Yll is symmetric) as
    //                   Yll * M' = Yhl'                 (3)
    // where, trivially, Yhl' = (Yh*Yl')' = Yl*Yh' since Yll is symmetric 
    // positive definite, we can use Cholesky decomposition of Yll to solve 
    // the linear system without explicitly inverting the matrix Yll, thus 
    // obtaining the MOST efficient of all the three strategies;
    // NOTE: the solution of (3) will be M transposed (!!!)
    #elif LIN_MAP_MODE==3
    // Yhl' <- Yl*Yh'
    gemm(REAL, 'N', 'T', s, s, n, Yl, s, Yh, s, Yhl, s);
    // Yll <- Yl*Yl'
    syrk(REAL, 'U', 'N', s, n, Yl, s,  Yll, s);
    // Yll <- Yll+lambda*I
    for(size_t d=0; d<s; d++) {
        Yll[d*(s+1)] += lambda;
    }
    // Cholesky factorization
    potrf('U', s, Yll, s);
    // Yhl <- (Yll^-1)*Yhl'  solve symmetric system, the solution (in this 
    //                       case M) is overwritten to the right hand side 
    //                       (Yhl); on exit Yhl will be equal to M transposed) 
    potrs('U', s, s, Yll, s, Yhl, s);
    // xh <- Yhl'*xl         the corrected patch is M*xl, being in this case M=Yhl'
    gemm(REAL, 'T', 'N', s, 1, s, Yhl, s, xl, s, xh, s);
    
    #else
    // wrong LIN_MAP_MODE parameter
    throw std::invalid_argument("linear_mapping : invalid LIN_MAP_MODE");
    #endif
    return xh;
}

#endif

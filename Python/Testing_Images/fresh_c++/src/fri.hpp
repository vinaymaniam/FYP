
#ifndef FRI_H
#define FRI_H

#include <algorithm>
#include <complex>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

#include "master.h"
#include "algebra.hpp"
#include "filters.hpp"

// unmodifiable parameter
#define S_M       "S_M"    // exponential moments
#define LSV_H     "LSV_H"  // left signular vectors of H
#define PINV_U    "PINV_U" // pseudo-inverse of U
#define EIG_Z     "EIG_Z"  // eigenvalues of Z

// modifiable parameter
#define POLYFIT    2 // integer {1->QR, 2->lin regression}
#define LEFT_EIG_H 1 // integer {1->SVD, 2->eigenproblem on H*H'}

class fri {
private:
    size_t N;     // length of 1-D signal
    size_t level; // level of decomposition
    real_t a;     // weight
    int filter;   // wavelet type {spln4 or bior4.4}
    
    real_t loc_shift;
    real_t T_s;
    std::vector<real_t> phi;
    
    std::vector<real_t> t;
    complex_t lambda;
    std::vector<real_t> c_m_n_mdf;
    
    #if PARALLEL==0 && VERBOSE>=3
    std::map<std::string, size_t> time_log;
    std::chrono::high_resolution_clock::time_point t_start, t_end;
    #endif
    
public:
    /*   default empty constructor */
    fri() {};
    
    /* constructor: computes internal vector t, constant lambda, and matrix  
       c_m_n_mdf from the length of 1-D signal N */
    fri(size_t _N, size_t _level, real_t _a, int _filter = SPLN4) {
        #if SAFE_CHECK
        if((_filter!=SPLN4 && _filter!=BIOR4_4) || _level>3 || _level==0) {
            throw std::invalid_argument("fri::fri : invalid parameters");
        }
        #endif
        // init map structures
        std::map<size_t, std::vector<std::vector<real_t> > > phi_map = PHI_MAP;
        std::map<size_t, std::vector<size_t> > loc_shift_map = LOC_SHIFT_MAP;
        
        // setting up private parameters
        N         = _N;
        level     = _level;
        a         = _a;
        filter    = _filter;
        
        #if PARALLEL==0 && VERBOSE>=3
        time_log = {{S_M,0}, {LSV_H,0}, {PINV_U,0}, {EIG_Z,0}};
        #endif
        
        // get wavelet filter (!!! shift = 0 !!!)
        phi = phi_map[filter][level-1];
        loc_shift = (real_t)loc_shift_map[filter][level-1];
        
        // start filter setup
        size_t P = N-1;         // order of exp repr kernel is P+1
        size_t L = 1<<level;    // 2^level
        real_t T = 1/(real_t)N; // time-scale factor
        
        T_s = T/(real_t)L; // time resolution
        int n1;
        if(N%2==0) {
            n1 = -(int)N/2;
        } else {
            n1 = -((int)N-1)/2;
        }
        // t definition
        real_t t0 = n1 * T;
        t = std::vector<real_t>(N*L);
        for(size_t i=0; i<t.size(); i++) {
            t[i] = t0 + T_s*(real_t)i;
        }
        // we assume to deal with sampling kernel phi that are able to
        // approximately reproduce exponential functions, formally
        //          sum_n c_m_n*phi(t-n) \approx e^(alpha_m*t)   (1)
        // where alpha_m = alpha_0+lambda*m, for m in [0,P]
        // in practice the best exponential reproducing kernels should 
        // reproduce exponentials with exponents of the form
        //               alpha_m = 1i*pi/(P+1)*(2m-P)            (2)
        // with m in [0,P] and the order of the kernel is such that N=P+1
        // however, in the approximate case, the denominator of (2) should 
        // be relaxed to allow the modulus of coefficients c_m_0 to be not  
        // exactly equal to 1 (i.e., the exact case), it can be proven that
        // the best values of the denominator L are such that P+1<=L<=4*(P+1)
        // the parameter a controls this relaxation and thus 1<=a<=4
        real_t L_den = a * (P+1);
        real_t alpha_0 = (real_t)(-M_PI * P / L_den);
        lambda = {0, (real_t)(2 * M_PI / L_den)};
        std::vector<complex_t> alpha_m;
        alpha_m.reserve(N);
        for(size_t i=0; i<N; i++) {
            alpha_m.push_back( {0, alpha_0 + lambda.imag() * (real_t)i} );
        }
        // generate B-Spline of zero-order (i.e., our assumed sampling kernel)
        std::vector<real_t> nphi(L+1, 1);
        nphi[0] = 0.5;
        nphi[nphi.size()-1] = 0.5;
        // in the approximate case, we are interested in the coefficients
        // c_m_n=c_m_0*e^(alpha_m*n) that best fit (1), which is solved by
        // rewriting (1) using the factorization of c_m_n as
        //   g_m_alpha(t) = c_m_0 * sum_n { e^(-alpha_m*(t-n))*phi(t-n) }
        // and then we can reduce the problem of finding c_m_n to the problem
        // of approximating the Fourier series expansion of g_m_alpha to 1
        // g_m_alpha(t) = sum_l { c_m_0*phi_hat(alpha_m+1i*2*pi*l)*e^(1i*pi*l*t) } \approx 1 (3)
        // thus, if phi only approximately reproduces exponential kernels
        // we have that the terms phi_hat(alpha_m+1i*2*pi*l) do not vanish,
        // however if phi_hat decays quickly, only few terms of (3) are 
        // needed to have an accurate bound for the error
        // there are different strategies to compute the coefficients c_m_0 
        // for an approximate reproduction: here we use the so-called 
        // constant least-squares strategy, where we naturally minimize the
        // reproduction error by choosing c_m_0=phi_hat(alpha_m)^-1, which
        // has the advantage of simplicity and also the advantage that
        // it only requires the knowledge of the Laplace transform of the 
        // sampling kernel phi at the locations alpha_m
        std::vector<real_t> newphi = T_s * conv1(phi, nphi, FULL);
        std::vector<complex_t> c_m_0;
        c_m_0.reserve(N);
        complex_t one = {1, 0};
        for(size_t i=0; i<N; i++) {
            complex_t laplace_coef = {0, 0};
            for(size_t j=0; j<newphi.size(); j++) {
                laplace_coef += std::exp( -((real_t)j/(real_t)L)*alpha_m[i] ) * newphi[j];
            }
            // invert the Laplace transform coefficients
            c_m_0.push_back( one / laplace_coef );
        }
        // calculate coefficients c_m_n_mdf (!!! shift = 0 !!!), i.e. the
        // matrix of coefficients that transforms the vector of samples
        // to the vector of moments s_m, by setting P=N-1, we make the
        // c_m_n_mdf matrix (P+1)-by-N a square matrix N-by-N
        // the coefficients in c_m_n_mdf depends on the sampling kernel
        // and note that in order to solve an FRI problem P+1>=2K, being K
        // the number of diracs that one needs to recover, thus making the
        // total degrees of freedom=2K (K locations and K amplitudes)
        
        // complex matrix N-by-N vectorized in column-major form using
        // FORTRAN convention, i.e. the ith complex element occupies the
        // positions 2*i and 2*i+1 in the vectorized matrix representing
        // real and imaginary parts, respectively
        c_m_n_mdf = std::vector<real_t>(N*N*COMPLEX);
        for(size_t j=0; j<N; j++) {
            real_t n = n1+(real_t)j;
            for(size_t i=0; i<N; i++) {
                // the entries to the c_m_n_mdf matrix are defined as:
                // phi_hat( alpha_m )^-1 * e^(alpha_m * n)
                complex_t coef = c_m_0[i] * std::exp(n*alpha_m[i]);
                c_m_n_mdf[COMPLEX*(j*N+i)]   = coef.real();
                c_m_n_mdf[COMPLEX*(j*N+i)+1] = coef.imag();
            }
        }
    }
    
    /* find optimal locations tk of discontinuities in the signal x using
       a threshold on the singluar values of the matrix H */
    std::vector<real_t> find_locations(const std::vector<real_t> &x, real_t thresh) {
        #if SAFE_CHECK
        // check if length of x is compatible with that of kernel
        if(c_m_n_mdf.empty() || N!=x.size()) {
            throw std::invalid_argument("fri::find_locations : invalid length of signal " +
                    std::to_string(x.size()) + " vs. kernel dimension " + std::to_string(N));
        }
        #endif
        
        // the coefficients of the exponential reproducing kernel c_m_n_mdf
        // can be adapted to particular classes of non-exponential functions
        // such as piece-wise constant signals (the ones assumed here, i.e.
        // any 1-D cross-section of natural images) by computing the finite
        // difference derivative of the signal samppled with phi, which are
        // equivalend to the measurements of the derivative of the signal
        // sampled with phi*beta (being beta a box function); note that
        // the derivative of the signal is a train of Diracs, thus we are
        // able to link the problem of piece-wise constant signal to the
        // known case of signals composed by a stream of Diracs
        
        #if PARALLEL==0 && VERBOSE>=3
        t_start = std::chrono::high_resolution_clock::now();
        #endif
        
        // create complex (FORTRAN) vector N-by-1 having real part equal to
        // a subtraction of one-point right shifted x and x (last elem = 0)
        // which represent the finite difference (derivative) of x
        std::vector<real_t> dx(x.size()*COMPLEX, 0);
        for(size_t i=0; i<N-1; i++) {
            dx[2*i] = x[i+1]-x[i];
        }
        // compute the vector of moments s_m by multiplying the coefficients
        // of exponential reproducing kernels with the input derivative dx
        std::vector<real_t> s_m(N*COMPLEX);
        gemm(COMPLEX, 'N', 'N', N, 1, N, c_m_n_mdf, N, dx, N, s_m, N);
        
        #if PARALLEL==0 && VERBOSE>=3
        t_end = std::chrono::high_resolution_clock::now();
        time_log.at(S_M) += MICROSEC(t_end-t_start);
        t_start = std::chrono::high_resolution_clock::now();
        #endif
        
        // the 2K unkowns (locations and amplitutes of the stream of Dirac)
        // can be recovered from s_m by the annihilation (Prony's) method 
        // by solving the system H*h_m=0, where H is a Toeplitz matrix 
        // m-by-n obtained from the moments s_m and h_m are the (unknown)
        // coefficients of the filter whose z-transform has roots uk that
        // can be easily related to the location tk as uk = e^(lambda*yk/T); 
        // the zeros uniquely defines uk if the locations tk are distinct;
        
        // create horizontally flipped Toeplitz matrix H from two subvectors
        // of the exponential moments s_m (H is roughly square):
        //         col=s_m[J:N-1]        row=reverse(s_m[0:J-1])
        // H is produced in colum-major vectorized form with real and imag
        // part of each element interleaved (suitable for LAPACK)
        size_t J = N/2 - (N%2==0);
        size_t m = s_m.size()/2-J;
        size_t n = s_m.size()/2-N+J+1;
        std::vector<real_t> H = std::vector<real_t>(m*n*COMPLEX);
        for(size_t j=0; j<n; j++) {
            size_t k = n-j-1;
            for(size_t i=0; i<m; i++) {
                // if i>j    select (i-j)-th element from col, i.e. the 
                //           (J+i-j)-th element from s_m
                // if i<=j   select (j-i)-th element from row, i.e. the 
                //           (J-j+i)-th element from s_m
                size_t idx_s_m = (J+i-j) * COMPLEX;
                size_t idx_H   = (k*m+i) * COMPLEX;
                H.at(idx_H)   = s_m.at(idx_s_m);
                H.at(idx_H+1) = s_m.at(idx_s_m+1);
            }
        }
        
        // explicitly, H*h_m=0 is s_m=sum_k{x_k*u_k^m}, where x_m and u_k
        // are the unknowns and we only know the values of s_m and that
        // the s_m's are in a power series form; the right null space of
        // H*h_m=0 is calculating the SVD of H and choosing the singular
        // vector corresponding to the zero singular value; however, when
        // the moments are noisy H*h_m=0 is not satisfied perfectly, H*h_m=0
        // is solved through least squares thus finding through the pseudo-
        // inverse the vector h_m that minimizes the l2-norm of the error, 
        // however when the moments s_m are not exact (noisy) a total-least-
        // squares (TLS) approach should be preferred because all complex 
        // moments are perturbed and error appear on both sides of the 
        // equations, i.e. the solution that minimizes the norm ||H*h_m||^2
        // under the constraint ||h_m||^2=1; this is solved by computing the
        // SVD of H=U*S*Vh, and then extracting the right singular vector
        // corresponding to the smallest singular values (the right-most 
        // column of V); however, we can do better, as we can solve the
        // structured TLS, which attempts to find the rank-K Toeplitz matrix
        // Hk that is closest (in some norm sense) to H using Cadzow's 
        // iterative algorithm; this is however not always converging to an 
        // optimum rank-K estimate; an alternative is directly estimate uk 
        // in a parametric fashion without the annihilating coefficients h_m
        // via the subspace estimator method which is based on the shift-
        // invariance property derived from the matrix pencil method which
        // finds directly uk as the eigenvalues of an operator that maps 
        // eithper U0 onto U1 or V0 onto V1, i.e., U1=U0*Z, where Z is
        // K-by-K matrix having uk as eigenvalues and U1 and U0 are matrices 
        // obtained after removing the first row (or column), and the last
        // row (or column) from U, respectively;
        
        // we are interested in the left eigenvectors of H, which we can 
        // compute either with an SVD on H, or as an eigenproblem on H*H'
        #if LEFT_EIG_H==1
        // calculate COMPLEX SVD of H in 'S' mode (i.e. only first MIN(m,n) 
        // columns of U and rows of Vt will be computed)
        std::vector<real_t> S(MIN(m,n));
        std::vector<real_t> U(m*m*COMPLEX);
        std::vector<real_t> Vt(n*n*COMPLEX);
        gesdd(COMPLEX, 'S', m, n, H, m, S, U, Vt);
        real_t sv_thresh = thresh;
        #elif LEFT_EIG_H==2 // <- default branch !!! 
        // calculate complex eigenproblem on H*H', because the eigenvectors
        // of H*H' are equal to the left eigenvectors U of the SVD
        // decomposition of H; also the non-zero signular values S of H are
        // equal to the square root of the eigenvalues of H*H'
        std::vector<real_t> U(m*m*COMPLEX);
        // multiply H by its transpose and store results to U
        syrk(COMPLEX, 'U', 'N', m, n, H, m, U, m);
        std::vector<real_t> S(m);
        // eigenproblem on symmetric matrix U=H*H', U is overwritten with
        // eigenvectors and S is overwritten with eigenvalues
        syevd(COMPLEX, 'V', 'U', m, U, m, S);
        // both eigenvectors and eigenvalues are stored in ascending order
        // thus we need to reverse it to be compliant with the ordering
        // of singular vectors/values obtained by SVD (which is expected
        // by the remainder of the algorithm)
        if(sum1(abs1(S))>EPSILON) {
            // do not reverse if all eigenvalues are zero (this can happen
            // if input signal x is constant, and thus H is all zero)
            for(size_t j=0; j<m/2; j++) {
              for(size_t i=0; i<COMPLEX*m; i+=2) {
                  std::swap( U[COMPLEX*j*m+i],   U[COMPLEX*(m-j-1)*m+i]   );
                  std::swap( U[COMPLEX*j*m+i+1], U[COMPLEX*(m-j-1)*m+i+1] );
              }
            }
            std::reverse(S.begin(), S.end());
        }
        real_t sv_thresh = SQUARE(thresh);
        #else
        throw std::invalid_argument("fri::find_locations : invalid LEFT_EIG_H");
        #endif
        // estimate K from singular values of H
        size_t K = 1;
        for(size_t i=1; i<MIN(m,n); i++) {
            K += (S[i]/S[0])>sv_thresh;
        }
        // work only on first K columns and the first n-p rows of U
        size_t p = 1;
        if(m<K || K==0) {
            #if WARNING
            std::cout << "Warning : fri::find_locations : changing " << 
                    "invalid K to the closest allowed value" << std::endl;
            #endif
            K = MIN(MAX(1,K),m);
        }
        // we can remove the first/last p rows of U, thus p<m
        if(p>=m) {
            throw std::out_of_range("fri::find_locations : invalid parameter p");
        }
        
        #if PARALLEL==0 && VERBOSE>=3
        t_end = std::chrono::high_resolution_clock::now();
        time_log.at(LSV_H) += MICROSEC(t_end-t_start);
        t_start = std::chrono::high_resolution_clock::now();
        #endif
        
        // U and Vt satisfies the shift-invariance property, thus U1=U0*Z,
        // where U1 is U after removing the first row, and U0 is U after
        // removing the last row, and Z is K-by-K diagonal matrix having
        // uk's as diagonal values; Z is found by solving U1=U0*Z via the
        // pseudoinverse of U0: Z=pinv(U0)*U1; note that U=[Uk Um-k+1], and 
        // the columns of Uk equivalently span the signal subspace, thus we
        // can focus only on the first K column of U, and find an operator
        // from the U0k and U1k truncated matrices of Uk: Z=pinv(Uk0)*Uk1;
        // note that in the presence of noise, the shif-invariance property
        // is not exaclty satisfied, however this has little effect on the
        // K dominant singular values and thus Z is still a good estimator
        
        // compute pseudo-inverse Upinv of U0=U[0:m-p-1,0:K-1] (m-p)-by-K, 
        // notice that U has still leading dimension equal to m
        // Upinv has size K-by-(m-p), i.e. equal to the transposed input
        std::vector<real_t> Upinv(K*(m-p)*COMPLEX);
        pinv(COMPLEX, m-p, K, U, m, Upinv);
        
        #if PARALLEL==0 && VERBOSE>=3
        t_end = std::chrono::high_resolution_clock::now();
        time_log.at(PINV_U) += MICROSEC(t_end-t_start);
        t_start = std::chrono::high_resolution_clock::now();
        #endif
        
        // product of pseudo-inverse Upinv K-by-(m-p) and U[p:m-1,0:K-1] (m-p)-by-K
        // thus the result Z=pinv(Uk0)*Uk1=Upinv*Uk1 has size K-by-K
        std::vector<real_t> Z(K*K*COMPLEX);
        // ignore first p rows of U from the multiplication using osb arg
        // if mode=COMPLEX we need to ignore two values (re,im) in each row
        gemm(COMPLEX, 'N', 'N', K, K, m-p, Upinv, K, U, m, Z, K, 0, p*COMPLEX);
        // calculate eigenvalues of the square matrix Z K-by-K
        std::vector<real_t> uk(K*COMPLEX);
        geev(COMPLEX, 'N', K, Z, K, uk);
        // apply power if p is different than one, the advantage of using 
        // values of p larger than one is that the separation among the 
        // estimated time delays is increased p times; this enhances the 
        // resolution capabilities of the original method, but also 
        // introduces ambiguities on the time delays of the uk
        if(p!=1) {
            uk = uk ^ (1.0f/(real_t)p);
        }
        
        #if PARALLEL==0 && VERBOSE>=3
        t_end = std::chrono::high_resolution_clock::now();
        time_log.at(EIG_Z) += MICROSEC(t_end-t_start);
        #endif
        
        // uk contains the K complex eigenvalues (the roots of the annihilation
        // filter) in FORTRAN format, i.e. the real and imaginari part of 
        // the k-th amplitute occupy positions 2*k and 2*k+1, respectively
        // from the amplitutes uk, we can retrieve the locations tk 
        std::vector<real_t> tk;
        tk.reserve(K);
        for(size_t i=0; i<K; i++) {
            // equivalent to real( log(uk) / (lambda*N) )
            // where eigenvalues and lambda are complex numbers, and lambda 
            // has real part equal to zero
            tk.push_back( std::atan2(uk[2*i+1],uk[2*i]) / (lambda.imag()*N) );
        }
        std::sort(tk.begin(), tk.end());
        
        return tk;
    }
    
    /* find location of discontinuities in x_low within a grid scaled by 
       factor equal to scale to build a signal of length n, the resulting
       discontinuities are given based on the threshold sv_thresh of the 
       singular values of the FRI model matrix */
    std::vector<real_t> discontinuities_location(const std::vector<real_t> &x_low, 
            size_t n, size_t scale, real_t sv_thresh) {
        #if SAFE_CHECK
        if(scale<=1) {
            throw std::invalid_argument("fri::discontinuities_location : invalid scale");
        }
        #endif
        
        // solve non-linear system to estimate optimal discontinuity positions
        // NOTE: when using float precision, d_pos might be different than
        //       the one obtained with double precision, as the numerical 
        //       errors previously accumulated (in particular due to the
        //       complex std::exp/log and xgemm operations) in number with 
        //       fractional part close to .5, might cause the round 
        //       operator to undershoot or overshoot the results in d_pos
        std::vector<real_t> d_pos = round1( 
                ( find_locations(x_low, sv_thresh) - (loc_shift*T_s+t[0]) ) * 
                ((real_t)scale / T_s) + (real_t)1.0 
        );
        
        // check that all discontinuities values represent a valid index in 
        // x_lin_up: note that d_pos is already sorted because tk is sorted
        // delete elements <=1 or >=n (x_lin_up.size())
        d_pos.erase( std::remove_if(d_pos.begin(), d_pos.end(), [&n](const real_t i) {
            return i<=1 || i>=n;
        }), d_pos.end() );
        // delete duplicates in d_pos
        d_pos.erase( std::unique(d_pos.begin(), d_pos.end()), d_pos.end() );
        // finally add 0 at the beginning
        d_pos.insert(d_pos.begin(), 0);
        // and add n at the and
        d_pos.push_back((real_t)n);
        // !!! sorting d_pos is not necessary because tk is already sorted
        return d_pos;
    }
    
    /* obtain upsampled signal of length equal to x_lin_up.size() where
       x_lin_up is a linear upsampling of the low resolution signal by 
       constructing a piecewise linear approximation using the locations
       in d_pos as the discontinuities in the signal and the values in 
       x_lin_up as the intensity used in the approximation; in other words
       for each pair of consecutive discontinuities in d_pos (d_i,d_{i+1})
       we extract the values v=x_lin_up(d_i:d_{i+1}) and we compute a 
       polyfit linear approximation on v, which will be the estimated
       upsampled signal */
    std::vector<real_t> upsample_line(const std::vector<real_t> &x_lin_up,
            const std::vector<real_t> &d_pos) {
        #if SAFE_CHECK
        if(d_pos.back()>x_lin_up.size()) {
            throw std::invalid_argument("fri::upsample_line : invalid locations");
        }
        #endif
        // fit polynomials onto x_lin_up between each pair of consecutive  
        // discontinuity positions in d_pos and store results in x_fri_up
        std::vector<real_t> x_fri_up(x_lin_up.size());
        for(size_t i=0; i<d_pos.size()-1; i++) {
            // starting x position (i-th discontinuity in d_pos)
            size_t x0 = (size_t)d_pos[i];
            // length of signal (range between i-th and i-th+1 d_pos
            size_t M = (size_t)d_pos[i+1]-x0;
            // for each interval between consecutive values in d_pos
            // provided that the interval contains more than one sample
            if(M>1) {
                // fit a polynomial of degree d to the range between each
                // consecutive pair of d_pos as x samples, and the values
                // in corresponding x positions in x_lin_up as y samples
                // the degree of polynomial is 1
                #if POLYFIT==1
                size_t d = 1;
                size_t D = d+1;
                // extract samples x_lin_up in range [x0:x0+M)
                std::vector<real_t> p(x_lin_up.begin()+x0, x_lin_up.begin()+x0+M);
                // construct Vandermonde matrix of the x positions, i.e., 
                // build model vectorized (column-major) matrix M-by-D with
                // each column i in {0, 1, ..., D} contains the (D-i)-th
                // power of the x samples, i.e., if x=[0 1 ... M] and D=d+1
                // [0^d 1^d ... M^d; 0^(d-1) ... M^(d-1); ... ; 0^0 ... M^0]
                // thus each column represent the contribution of a particular
                // power, being the last one the contribution of the constant
                // this is actually implemeted by initializing a matrix of
                // ones, and then iteratively multiplying (i-1)-th column
                // with the i-th column, with i ranging from 0 to D-2
                std::vector<real_t> V(M*D, 1);
                for(size_t i=0; i<D-1; i++) {
                    for(size_t j=0; j<M; j++) {
                        V.at((d-i-1)*M+j) = (x0+j)*V.at((d-i)*M+j);
                    }
                }
                // fitting data to a polynomial of degree d is a linear 
                // least square problem, which we solve here by first 
                // calculating the QR factorization of the model matrix V
                // which is actually implemeted by LAPACK functions geqrf
                // which returns the D-by-1 scalar factors tau; V on exit  
                // is overwritten in the upper-triangular part by R, and in 
                // the lower-triangular part by the reflectors that, combined 
                // with the factors tau, produce the matrix Q
                std::vector<real_t> tau(D);
                geqrf(M, D, V, M, tau);
                // here we compute the product Q'*p, being Q represented
                // by the reflectors contained in V and their factors tau;
                // on exit, the result of size D-by-1 in overwritten in p
                ormqr('L', 'T', M, 1, D, V, M, tau, p, M);
                // finally we solve the linear system R*x=p, where R is in
                // the upper triangular part of V; on exit, p is the result
                trtrs('U', 'N', 'N', D, 1, V, M, p, M);
                #elif POLYFIT==2
                // instead of solving least-squares problem via QR, let us
                // simply solve the fitting via linear regression, i.e. his
                // is achieved by calculating the derivative with respect to
                // the polynomial (line) parameters and setting these to zero
                // average of x and y data
                real_t x_bar = x0+((real_t)M-1)/2;
                real_t y_bar = 0;
                for(size_t i=0; i<M; i++) {
                    y_bar += x_lin_up[x0+i];
                }
                y_bar /= M;
                // derivatives 
                real_t Sxy = 0;
                real_t Sx = 0;
                for(size_t i=0; i<M; i++) {
                    real_t dx = x0+i-x_bar;
                    real_t dy = x_lin_up[x0+i]-y_bar;
                    Sx  += SQUARE(dx);
                    Sxy += dx*dy;
                }
                // regression coefficients p = {slope, intercept}
                std::vector<real_t>p(2);
                p[0] = Sxy/Sx;
                p[1] = y_bar-p[0]*x_bar;
                #else
                // wrong POLYFIT parameter
                throw std::invalid_argument("fri::upsample_line : invalid POLYFIT");
                #endif
                // polynomial coefficients are the first N elements in p
                for(size_t i=0; i<M; i++) {
                    x_fri_up.at(x0+i) = p[0] * (x0+i) + p[1];
                }
            } else {
                // if interval has length 1
                x_fri_up.at(x0) = x_lin_up.at(x0);
            }
        }
        
        // linfit is upplied to scaled version of signal, so remember to 
        // compute the  wavelet decomposition on the fitted line x_fri_up 
        // to get the final result of the expected size
        return x_fri_up;
    }
    
    #if PARALLEL==0 && VERBOSE>=3
    const std::map<std::string, size_t>& get_time_log() const {
        return time_log;
    }
    #endif
    
};

#endif

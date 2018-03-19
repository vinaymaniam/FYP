
#ifndef FOURIER_H
#define FOURIER_H

#include <complex>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "master.h"
#include "utils.hpp"
#include "algebra.hpp"
#include "matrix.hpp"
#include "volume.hpp"

// unmodifiable constants
#define FORWARD 0
#define INVERSE 1

class fourier {
private:
    std::unordered_map<size_t, std::vector<real_t> > dft_matrices;
    typedef std::vector<complex_t> (fourier::*dft_1d_op)(const std::vector<complex_t> &in, int dir);
    
public:
    /* default empty constructor */
    fourier() {}
    
    /* DFT on 1-D complex input using complex matrix multiplication */
    std::vector<complex_t> dft_1d_transform(const std::vector<complex_t> &in, int dir) {
        #if SAFE_CHECK
        if(dir!=FORWARD && dir!=INVERSE) {
            throw std::invalid_argument("fourier::dft_1d_transform : invalid operation");
        }
        #endif
        size_t N = in.size();
        char trans_dft_mtx = (dir==FORWARD)? 'N' : 'C';
        std::vector<real_t> out(N*2);
        gemm(COMPLEX, 'N', trans_dft_mtx, 1, N, N, c2fort(in), 1, fftmtx(N), N, out, 1);
        return (dir==FORWARD)? fort2c(out) : fort2c(out / (real_t)N);
    }
    
    /* 1-D fftshift operation */
    std::vector<complex_t> dft_1d_shift(const std::vector<complex_t> &in, int dir) {
        #if SAFE_CHECK
        if(dir!=FORWARD && dir!=INVERSE) {
            throw std::invalid_argument("fourier::dft_1d_shift : invalid operation");
        }
        #endif
        size_t N = in.size();
        size_t h = std::ceil((real_t)N/2) - (dir==INVERSE && N%2==1);
        std::vector<complex_t> out;
        out.reserve(N);
        for(size_t i=0; i<N; i++) {
            out.push_back( in[(h+i)%N] );
        }
        return out;
    }
    
    /* apply separable 1-D operation f_id on 2-D input */
    matrix<complex_t> dft_2d_op_dispatcher(dft_1d_op f_1d, const matrix<complex_t> &in, int dir) {
        #if SAFE_CHECK
        if(dir!=FORWARD && dir!=INVERSE) {
            throw std::invalid_argument("fourier::dft_2d_op_dispatcher : invalid operation");
        }
        #endif
        matrix<complex_t> out = matrix<complex_t>(in.size());
        // apply 1-D operation on each row of input
        for(size_t i=0; i<in.getM(); i++) {
            out.set_row( (this->*f_1d)(in.get_row(i), dir), i );
        }
        // apply 1-D operation on each column of partially processed output
        for(size_t i=0; i<in.getN(); i++) {
            out.set_col( (this->*f_1d)(out.get_col(i), dir), i );
        }
        return out;
    }
    
    /* apply separable 1-D operation f_id on 3-D input */
    volume<complex_t> dft_3d_op_dispatcher(dft_1d_op f_1d, const volume<complex_t> &in, int dir) {
        #if SAFE_CHECK
        if(dir!=FORWARD && dir!=INVERSE) {
            throw std::invalid_argument("fourier::dft_3d_op_dispatcher : invalid operation");
        }
        #endif
        volume<complex_t> out = volume<complex_t>(in.size());
        // for each plane of input
        for(size_t k=0; k<in.getP(); k++) {
            // apply 1-D operation on each row of input
            for(size_t i=0; i<in.getM(); i++) {
                out.set_row( (this->*f_1d)(in.get_row(i,k), dir), i, k );
            }
            // apply 1-D operation on each column of partially processed output
            for(size_t i=0; i<in.getN(); i++) {
                out.set_col( (this->*f_1d)(out.get_col(i, k), dir), i, k );
            }
        }
        // apply 1-D operation on each cross-section
        for(size_t j=0; j<out.getN(); j++) {
            for(size_t i=0; i<out.getM(); i++) {
                out.set_xsec( (this->*f_1d)(out.get_xsec(i,j),dir), i, j );
            }
        }
        return out;
    }
    
    /* 1-D DFT */
    std::vector<complex_t> fft(const std::vector<complex_t> &in) {
        return dft_1d_transform(in, FORWARD);
    }
    
    /* inverse 1-D DFT */
    std::vector<complex_t> ifft(const std::vector<complex_t> &in) {
        return dft_1d_transform(in, INVERSE);
    }
    
    /* 2-D DFT */
    matrix<complex_t> fft2(const matrix<complex_t> &in) {
        return dft_2d_op_dispatcher(&fourier::dft_1d_transform, in, FORWARD);
    }
    
    /* inverse 2-D DFT */
    matrix<complex_t> ifft2(const matrix<complex_t> &in) {
        return dft_2d_op_dispatcher(&fourier::dft_1d_transform, in, INVERSE);
    }
    
    /* 3-D DFT */
    volume<complex_t> fft3(const volume<complex_t> &in) {
        return dft_3d_op_dispatcher(&fourier::dft_1d_transform, in, FORWARD);
    }
    
    /* inverse 3-D DFT */
    volume<complex_t> ifft3(const volume<complex_t> &in) {
        return dft_3d_op_dispatcher(&fourier::dft_1d_transform, in, INVERSE);
    }
    
    /* 1-D fftshift */
    std::vector<complex_t> fftshift(const std::vector<complex_t> &in) {
        return dft_1d_shift(in, FORWARD);
    }
    
    /* inverse 1-D fftshift */
    std::vector<complex_t> ifftshift(const std::vector<complex_t> &in) {
        return dft_1d_shift(in, INVERSE);
    }
    
    /* 2-D fftshift */
    matrix<complex_t> fftshift2(const matrix<complex_t> &in) {
        return dft_2d_op_dispatcher(&fourier::dft_1d_shift, in, FORWARD);
    }
    
    /* inverse 2-D fftshift */
    matrix<complex_t> ifftshift2(const matrix<complex_t> &in) {
        return dft_2d_op_dispatcher(&fourier::dft_1d_shift, in, INVERSE);
    }
    
    /* 3-D fftshift */
    volume<complex_t> fftshift3(const volume<complex_t> &in) {
        return dft_3d_op_dispatcher(&fourier::dft_1d_shift, in, FORWARD);
    }
    
    /* inverse 3-D fftshift */
    volume<complex_t> ifftshift3(const volume<complex_t> &in) {
        return dft_3d_op_dispatcher(&fourier::dft_1d_shift, in, INVERSE);
    }
    
    /* get DFT matrix of size N-by-N */
    std::vector<real_t> fftmtx(size_t N) {
        auto dft_matrix_iter = dft_matrices.find(N);
        if(dft_matrix_iter==dft_matrices.end()) {
            // N-by-N DFT matrix not found, so let us create a new one
            std::complex<real_t> omega = {std::cos(-2*M_PI/N), std::sin(-2*M_PI/N)};
            std::vector<real_t> dft_matrix;
            dft_matrix.reserve(N*N*2);
            for(size_t j=0; j<N; j++) {
                for(size_t i=0; i<N; i++) {
                    complex_t omega_ij = std::pow(omega,i*j);
                    dft_matrix.push_back(omega_ij.real());
                    dft_matrix.push_back(omega_ij.imag());
                }
            }
            dft_matrices.insert(std::make_pair( N, dft_matrix ));
            return dft_matrix;
        } else {
            return dft_matrix_iter->second;
        }
    }
    
};

#endif

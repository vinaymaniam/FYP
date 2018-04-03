
#ifndef IMAGE_H
#define IMAGE_H

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "master.h"
#include "algebra.hpp"
#include "utils.hpp"
#include "matrix.hpp"
#include "numeric.hpp"

// unmodifiable parameters
#define NO_ROUND   0
#define ROUND      1
#define BILINEAR   0
#define BICUBIC    1
#define DIR_X      0
#define DIR_Y      1

typedef real_t (*kernel) (real_t);

template<class T> using image_t = std::vector<matrix<T> >;

/* transform color space from RGB to YCbCr, using transformation matrix from 
   equations (9.6) in Poynton's, "Introduction to Digital Video" (p. 176)
   input and output are vector of matrix object, each representing a channel 
   in the corresponding image, i.e. {R, G, B} and {Y, Cb, Cr} */
template<class T>
image_t<T> rgb2ycbcr(const image_t<T> &rgb, int mode = NO_ROUND) {
    #if SAFE_CHECK
    if(rgb.size()!=3) {
        throw std::invalid_argument("rgb2ycbcr : invalid number of channels in input");
    }
    #endif
    // column-major vectorized transform matrix T (R, G, B are assumed to be in [0,255])
    std::vector<real_t> TM = {65.481f, -37.797f, 112.0f, 128.553f, -74.203f, -93.786f, 24.966f, 112.0f, -18.2140f};
    TM = TM / (real_t)255.0;
    // offset for each channel
    std::vector<real_t> offset = {16, 128, 128};
    // initialize output image
    std::vector<matrix<T> > yuv;
    yuv.reserve(rgb.size());
    for(size_t c=0; c<rgb.size(); c++) {
        // conversion: YCbCr = RGB * TM + offset
        yuv.push_back( rgb[0] * TM[c] + rgb[1] * TM[3+c] + rgb[2] * TM[6+c] + offset[c] );
        // round values if mode is set to ROUND
        if(mode==ROUND) {
            yuv[c].round();
        }
    }
    return yuv;
}

/* transform color space from YCbCr to RGB, using transformation matrix from 
   equations (9.6) in Poynton's, "Introduction to Digital Video" (p. 176)
   input and output are vector of matrix object, each representing a channel 
   in the corresponding image, i.e. {Y, Cb, Cr} and {R, G, B} */
template<class T>
image_t<T> ycbcr2rgb(const image_t<T> &yuv, int mode = NO_ROUND) {
    #if SAFE_CHECK
    if(yuv.size()!=3) {
        throw std::invalid_argument("ycbcr2rgb : invalid number of channels in input");
    }
    #endif
    // column-major vectorized inverse transform matrix T (channels assumed to be in [0,255])
    std::vector<real_t> TMinv = {
        0.0045662100456621f, 0.0045662100456621f, 0.0045662100456621f, 
        0.0000000011808800f, -0.0015363236860449f, 0.0079107162335547f,
        0.0062589289699439f, -0.0031881109496557f, 0.0000000119774970f
    };
    // offset for each channel TMinv*[16 128 128]'
    std::vector<real_t> offset_inv = {0.8742024200360562f, -0.5316682726390843f, 1.0856325717452218f};
    // scale TMinv and offset_inv
    TMinv = TMinv * (real_t)255.0;
    offset_inv = offset_inv * (real_t)255.0;
    // initialize output image
    image_t<T> rgb;
    rgb.reserve(yuv.size());
    for(size_t c=0; c<yuv.size(); c++) {
        // conversion: RGB = Tinv * YCbCr - offset_inv
        rgb.push_back( yuv[0] * TMinv[c] + yuv[1] * TMinv[3+c] + yuv[2] * TMinv[6+c] - offset_inv[c] );
        // round values if mode is set to ROUND
        if(mode==ROUND) {
            rgb[c].round();
        }
    }
    return rgb;
}

/* calculate gradient magnitude and direction */
template<class T>
void gradient(const matrix<T> &x, matrix<T> &Gmag, matrix<T> &Gdir) {
    #if SAFE_CHECK
    if(Gmag.getM()<x.getM() || Gmag.getN()<x.getN() || Gdir.getM()<x.getM() || Gdir.getN()<x.getN()) {
        throw std::invalid_argument("gradient : invalid output dimensions");
    }
    #endif
    // sobel kernel h1*h2
    std::vector<T> h1 = {1, 0, -1};
    std::vector<T> h2 = {1, 2, 1};
    // calculate horizontal and vertical gradients
    matrix<T> Gx = gradientxy(x, h1, h2, DIR_X);
    matrix<T> Gy = gradientxy(x, h1, h2, DIR_Y);
    // calculate magnitude and direction of gradient
    for(size_t i=0; i<x.getM(); i++) {
        for(size_t j=0; j<x.getN(); j++) {
            // hypot
            Gmag(i,j) = std::sqrt( SQUARE(std::abs(Gx(i,j))) + SQUARE(std::abs(Gy(i,j))) );
            // convert radian to degree
            Gdir(i,j) = std::atan2(-Gy(i,j),Gx(i,j))*180/M_PI;
        }
    }
}

template<class T>
matrix<T> gradientxy(const matrix<T> &x, const std::vector<T> &h1, const std::vector<T> &h2, int dir = DIR_X) {
    // create function pointers (default DIR_X)
    typedef std::vector<T> (matrix<T>::*getter) (size_t) const;
    typedef void (matrix<T>::*setter) (const std::vector<T>&, size_t);
    getter get1 = &matrix<T>::get_row;
    getter get2 = &matrix<T>::get_col;
    setter set1 = &matrix<T>::set_row;
    setter set2 = &matrix<T>::set_col;
    size_t dim1 = x.getM();
    size_t dim2 = x.getN();
    int xpad = 2;
    int ypad = 0;
    if(dir==DIR_Y) {
        std::swap(get1,get2);
        std::swap(set1,set2);
        std::swap(dim1,dim2);
        std::swap(xpad,ypad);
    } 
    // calculate gradient along one direction by 2D convolution with h1*h2
    matrix<T> tmp = matrix<T>(x.getM()+xpad,x.getN()+ypad);
    for(size_t i=0; i<dim1+2; i++) {
        // padding by replication along the vertical dir for pad=1
        //  (i<=pad)? 0 : (i>=dim)? dim-1 : i-pad;
        size_t k = (i<=1)? 0 : (i>=dim1)? dim1-1 : i-1;
        std::vector<T> xi = (x.*get1)(k);
        // padding by replication along the horizontal dir for pad=1
        xi.insert(xi.begin(),xi[0]);
        xi.push_back(xi.back());
        std::vector<T> y = conv1(xi,h1,VALID);
        (tmp.*set1)(y,i);
    }
    matrix<T> G = matrix<T>(x.getM(),x.getN());
    for(size_t i=0; i<dim2; i++) {
        std::vector<T> xi = (tmp.*get2)(i);
        std::vector<T> y = conv1(xi,h2,VALID);
        (G.*set2)(y,i);
    }
    return G;
}

/* kernel for cubic convolutional resampling */
inline real_t cubic_kernel(real_t x) {
    real_t absx = std::abs(x);
    real_t absx2 = absx*absx;
    real_t absx3 = absx2*absx;
    
    return (real_t)( (1.5*absx3 - 2.5*absx2 + 1) * (absx<=1) +
            (-0.5*absx3 + 2.5*absx2 - 4*absx + 2) * ((1<absx) & (absx<=2)) );
}

/* kernel for linear convolutional resampling */
inline real_t triangle_kernel(real_t x) {
    return (x+1) * ((-1<=x) & (x<0)) + 
            (1-x) * ((0<=x) & (x<=1));
}

/* computes and normilises convolution weights for the specified params */
inline std::vector<real_t> get_conv_weights(real_t u, real_t l, 
        real_t kernel_scale, const kernel &h, size_t conv_length) {
    std::vector<real_t> w(conv_length);
    real_t sum_w = 0;
    for(size_t k=0; k<conv_length; k++) {
        w[k] = kernel_scale*h( kernel_scale*(u-l-k) );
        sum_w += w[k];
    }
    return w / sum_w;
}

/* resize matrix x from M-by-N to ceil(M*scale)-by-ceil(N*scale) */
template<class T>
matrix<T> resize(const matrix<T> &x, real_t scale, int mode = BICUBIC) {
    return resize(x, scale, scale, mode);
}

/* resize matrix x from M-by-N to ceil(M*scaleM)-by-ceil(N*scaleN) */
template<class T>
matrix<T> resize(const matrix<T> &x, real_t scaleM, real_t scaleN, int mode = BICUBIC) {
    size_t M  = x.getM();
    size_t N  = x.getN();
    size_t Ms = (size_t)std::ceil(M*scaleM);
    size_t Ns = (size_t)std::ceil(N*scaleN);
    
    // setting up kernel (default is BICUBIC)
    size_t kernel_width = 4;
    kernel h = &(cubic_kernel);
    if(mode==BILINEAR) {
        kernel_width = 2;
        h = &(triangle_kernel);
    }
    
    // scale the kernel for resizing columns
    real_t kernel_scale = (scaleM<1)? scaleM : 1;
    real_t scaled_kernel_width = (real_t)kernel_width/kernel_scale;
    size_t conv_length = (size_t)std::ceil(scaled_kernel_width) + 2;
    
    // resizing columns
    matrix<T> tmp = matrix<T>(Ms,N);
    for(size_t i=0; i<Ms; i++) {
        real_t u = (i+1)/scaleM + 0.5f*(1-1/scaleM);
        real_t l = std::floor(u - scaled_kernel_width/2);
        // first the weights for ith row are computed
        std::vector<real_t> w = get_conv_weights(u, l, kernel_scale, 
                h, conv_length);
        // then we convolve w with ith row
        for(size_t k=0; k<conv_length; k++) {
            size_t p = (size_t)MIN( MAX(0, l+k-1), M-1 );
            for(size_t j=0; j<N; j++) {
                tmp(i,j) = tmp(i,j) + w[k]*x(p,j);
            }
        }
    }
    
    // scale the kernel for resizing rows
    kernel_scale = (scaleN<1)? scaleN : 1;
    scaled_kernel_width = (real_t)kernel_width/kernel_scale;
    conv_length = (size_t)std::ceil(scaled_kernel_width) + 2;
    
    // resizing rows
    matrix<T> y = matrix<T>(Ms,Ns);
    for(size_t j=0; j<Ns; j++) {
        real_t u = (j+1)/scaleN + 0.5f*(1-1/scaleN);
        real_t l = std::floor(u - scaled_kernel_width/2);
        std::vector<real_t> w = get_conv_weights(u, l, kernel_scale, 
                h, conv_length);
        for(size_t k=0; k<conv_length; k++) {
            size_t p = (size_t)MIN( MAX(0, l+k-1), N-1 );
            for(size_t i=0; i<Ms; i++) {
                y(i,j) = y(i,j) + w[k]*tmp(i,p);
            }
        }
    }
    return y;
}

/* return a 2-D gaussian window of size N-by-N and standard deviation sigma */
std::vector<real_t> gausswin2d(real_t N, real_t sigma) {
    std::vector<real_t> w1d(N, 0);
    std::vector<real_t> w2d(N*N, 0);
    // build 1-D gaussian function
    real_t Nhalf = (N-1)/2;
    for(real_t i=0; i<N; i++) {
        w1d.at(i) = std::exp( -0.5 * SQUARE( sigma*(i-Nhalf)/Nhalf ) );
    }
    // build 2-D window by vector multiplication
    gemm(REAL, 'N', 'T', N, N, 1, w1d, N, w1d, N, w2d, N);
    return w2d;
} 

#endif

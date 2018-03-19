
#ifndef NUMERIC_H
#define NUMERIC_H

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

#include "master.h"
#include "algebra.hpp"

// unmodifiable parameters
#define FULL      0
#define SAME      1
#define VALID     2

/* compute squared l2-norm of the difference between x1 and x2; if the sizes
   of the inputs are different, the result will be the difference of the 
   subvectors indexed from 0 and min(length(x1),length(x2)) */
template<class T>
inline T squared_euclidean_distance(const std::vector<T> &x1, const std::vector<T> &x2) {
    #if SAFE_CHECK
    if(x1.size()==0 || x2.size()==0) {
        throw std::invalid_argument("squared_euclidean_distance : empty input");
    }
    #endif
    T l2norm = 0;
    size_t p;
    size_t s = MIN(x1.size(),x2.size());
    for(p=0; p<=s-4; p+=4) {
        T err0 = x1[p]-x2[p];
        T err1 = x1[p+1]-x2[p+1];
        T err2 = x1[p+2]-x2[p+2];
        T err3 = x1[p+3]-x2[p+3];
        l2norm += SQUARE(err0) + SQUARE(err1) + SQUARE(err2) + SQUARE(err3);
    }
    for( ; p<s; p++) {
        real_t err = x1[p]-x2[p];
        l2norm += SQUARE(err);
    }
    return l2norm;
}

/* find prime factors of the input */
template<class T>
std::vector<size_t> factorize(T in) {
    size_t n = size_t(in);
    std::vector<size_t> f;
    f.push_back(1);
    size_t z = 2;
    while(z<=n) {
        if(n%z==0) {
            f.push_back( z );
            n /= z;
        } else {
            z++;
        }
    }
    return f;
}

/* 1D convolution (full is default) */
template<class T>
std::vector<T> conv1(const std::vector<T> &x, const std::vector<T> &h, int mode = FULL) {
    size_t xsz = x.size();
    size_t hsz = h.size();
    size_t ysz = xsz+hsz-1;
    #if SAFE_CHECK
    if(hsz>xsz) {
        throw std::invalid_argument("conv1 : invalid kernel size");
    }
    #endif
    // compute convolution
    std::vector<T> out(ysz, 0);
    for(size_t i=0; i<out.size(); i++) {
        size_t jmin = (i<hsz-1)?   0 : i-hsz+1;
        size_t jmax = (i<xsz-1)? i+1 : xsz;
        for(size_t j=jmin; j<jmax; j++) {
            out[i] += x[j]*h[i-j];
        }
    }
    // extract subvector depending on mode
    std::vector<T> y;
    if(mode==SAME) {
        y.reserve(xsz);
        size_t p = hsz/2;
        y.insert(y.begin(), out.begin()+p, out.end()-p+(hsz%2==0));
    } else if(mode==VALID) {
        y.reserve(out.size()-2*hsz+2);
        y.insert(y.begin(), out.begin()+hsz-1, out.end()-hsz+1);
    } else {
        y = out;
    }
    return y;
}

/* 1-D piecewise cubic Hermite spline interpolation, tangents are 
   calculated as Monotone with Lam harmonic mean */
template<class T>
std::vector<T> interp1_pchip(const std::vector<T> &x, const std::vector<T> &y, const std::vector<T> &xi) {
    #if SAFE_CHECK
    if(x.size()!=y.size()) {
        throw std::invalid_argument("interp1_pchip : invalid dimensions");
    }
    #endif
    size_t n = x.size();
    std::vector<T> yp(n,0);
    std::vector<T> h(n-1,0);
    std::vector<T> d(n-1,0);
    // slopes yp(i) at interior points:
    //  average of d(i-1) and d(i) when they have the same sign.
    //  0 when d(i-1) and d(i) have opposites signs or either is zero.
    for(size_t i=0; i<n-1; i++) {
        h[i] = x[i+1] - x[i];
        d[i] = (y[i+1] - y[i])/h[i];
    }
    if(n==2) {
        yp.assign(n,d[0]);
    } else {
        for(size_t i=1; i<n-1; i++) {
            if(d[i-1]*d[i]>0) {
                yp[i] = 2*d[i-1]*d[i]/(d[i-1] + d[i]);
            } else {
                yp[i] = 0;
            }
        }
        if(d[0]*(2*d[0]-yp[1])>0) {
            yp[0] = 2*d[0] - yp[1];
        } else {
            yp[0] = 0;
        }
        if(d[n-2]*(2*d[n-2]-yp[n-2])>0) {
            yp[n-1] = 2*d[n-2] - yp[n-2];
        } else {
            yp[n-1] = 0;
        }
        // slopes at end points:
        //  yp(0) and yp(n-1) as non-centered, shape-preserving three-point formulae.
        yp[0] = ((2*h[0]+h[1])*d[0] - h[0]*d[1])/(h[0]+h[1]);
        if(yp[0]*d[0]<0) {
            yp[0] = 0;
        } else if((d[0]*d[1]<0) && (std::abs(yp[0])>std::abs(3*d[0]))) {
            yp[0] = 3*d[0];
        }
        yp[n-1] = ((2*h[n-2]+h[n-3])*d[n-2] - h[n-2]*d[n-3])/(h[n-2]+h[n-3]);
        if(yp[n-1]*d[n-2]<0) {
            yp[n-1] = 0;
        } else if((d[n-2]*d[n-3]<0) && (std::abs(yp[n-1])>std::abs(3*d[n-2]))) {
            yp[n-1] = 3*d[n-2];
        }
    }
    std::vector<T> yi(xi.size(),0);
    for(size_t i=0; i<xi.size(); i++) {
        size_t klo = 0;
        size_t khi = n-1;
        while(khi-klo>1) {
            size_t k = (size_t)FIX((khi+klo)/2.0);
            if(x[k]>xi[i]) {
                khi = k;
            } else {
                klo = k;
            }
        }
        if(khi>=n || klo>=n) {
            throw std::logic_error("interp1_pchip : invalid pivot samples");
        }
        real_t dx = x[khi] - x[klo];
        if(dx==0.0) {
            throw std::logic_error("interp1_pchip : sample points must be distinct");
        }
        T a = yp[klo];
        T b = yp[khi];
        if(y[khi]==y[klo]) {
            a = 0;
            b = 0;
        } else {
            T alpha = a / ( (y[khi] - y[klo])/dx );
            T beta  = b / ( (y[khi] - y[klo])/dx );
            if(alpha<0 || beta<0) {
                a = 0;
                b = 0;
            } else {
                if((alpha*alpha + beta*beta)>9) {
                    T tau = 3/sqrt(alpha*alpha + beta*beta);
                    a = tau*alpha*(y[khi] - y[klo])/dx;
                    b = tau*beta*(y[khi] - y[klo])/dx;
                }
            }
        }
        // evaluate cubic Hermite polynomial
        T t   = (xi[i] - x[klo])/dx;
        T t2  = t*t;
        T t3  = t2*t;
        T h00 = 2*t3 - 3*t2 + 1;
        T h10 = t3 - 2*t2 + t;
        T h01 = -2*t3 + 3*t2;
        T h11 = t3 - t2;
        // interpolated value
        yi[i] = h00*y[klo] + h10*dx*a + h01*y[khi] + h11*dx*b;
    }
    return yi;
}

/* equivalent to 'v5cubic' interpolation in Matlab */
template<class T>
std::vector<T> interp1_cubic(const std::vector<T> &x, const std::vector<T> &y, 
        const std::vector<T> &xi, T extrapval = NAN) {
    #if SAFE_CHECK
    if(x.size()!=y.size()) {
        throw std::invalid_argument("interp1_pchip : invalid dimensions");
    }
    #endif
    size_t n = x.size();
    std::vector<T> yi(xi.size(), 0);
    for(size_t i=0; i<xi.size(); i++) {
        if(xi[i]<x[0] || xi[i]>x[n-1]) {
            yi[i] = extrapval;
            continue;
        }
        size_t klo = 0;
        size_t khi = n-1;
        while(khi>klo+1) {
            size_t k = (size_t)FIX((khi+klo)/2.0);
            if(x[k]>xi[i]) {
                khi = k;
            } else {
                klo = k;
            }
        }
        if(khi>=n || klo>=n) {
            throw std::logic_error("interp1_cubic : invalid pivot samples");
        }
        T dx = x[khi] - x[klo];
        if(dx==0.0) {
            throw std::logic_error("interp1_cubic : sample points must be distinct");
        }
        size_t km1 = klo - 1;
        size_t kp1 = khi + 1;
        T a = (klo<1)?   3*y[klo] - 3*y[khi] + y[khi+1] : y[km1];
        T d = (kp1>n-1)? 3*y[khi] - 3*y[klo] + y[klo-1] : y[kp1];
        T b = y[klo];
        T c = y[khi];
        // evaluate cubic polynomial
        T t = (xi[i] - x[klo])/dx;
        T t2 = t*t;
        T t3 = t2*t;
        T c00 = (-t3 + 2*t2 - t)/2.0f;
        T c10 = (3*t3 - 5*t2 + 2)/2.0f;
        T c20 = (-3*t3 + 4*t2 + t)/2.0f;
        T c30 = (t3 - t2)/2.0f;
        // interpolated value
        yi[i] = a*c00 + b*c10 + c*c20 + d*c30;
    }
    return yi;
}

/* linear interpolation */
template<class T>
std::vector<T> interp1_linear(const std::vector<T> &x, const std::vector<T> &y, 
        const std::vector<T> &xi, T extrapval = NAN) {
    #if SAFE_CHECK
    if(x.size()!=y.size()) {
        throw std::invalid_argument("interp1_pchip : invalid dimensions");
    }
    #endif
    std::vector<size_t> xi_pos(xi.size(),0);
    size_t k = 1;
    // find k of the closest element x[k] smaller than xi[i] for each i
    for(size_t i=0; i<xi.size(); i++) {
        while(x[k]<xi[i]) {
            k++;
        }
        xi_pos[i] = k-1;
    }
    std::vector<T> yi(xi.size(),0);
    for(size_t i=0; i<xi.size(); i++) {
        if(xi[i]<x[0] || xi[i]>x[x.size()-1]) {
           yi[i] = extrapval;
        }
        T dxi = xi[i]-x[xi_pos[i]];
        T dx  = x[xi_pos[i]+1]-x[xi_pos[i]];
        T t  = dxi/dx;
        // interpolated value
        yi[i] = y[xi_pos[i]] + t*(y[xi_pos[i]+1]-y[xi_pos[i]]);
    }
    return yi;
}

#endif


#ifndef WAVELET_H
#define WAVELET_H

#include <stdexcept>
#include <vector>

#include "master.h"
#include "numeric.hpp"
#include "matrix.hpp"
#include "filters.hpp"

// unmodifiable constants
#define LL_LH      0
#define HL_HH      1
#define DIAGONAL   0
#define VERTICAL   1
#define HORIZONTAL 2

template<class T> 
        class wavelet {
private:
    std::vector<T> Lo_D; // decomposition low-pass filter
    std::vector<T> Hi_D; // decomposition high-pass filter 
    std::vector<T> Lo_R; // reconstruction low-pass filter 
    std::vector<T> Hi_R; // reconstruction high-pass filter
    
    size_t length_filter;
    
public:
    /* constructor passing Lo and Hi dec and rec filters */
    wavelet(const std::vector<T> &_Lo_D, const std::vector<T> &_Hi_D, 
            const std::vector<T> &_Lo_R, const std::vector<T> &_Hi_R) {
        Lo_D = _Lo_D;
        Hi_D = _Hi_D;
        Lo_R = _Lo_R;
        Hi_R = _Hi_R;
        length_filter = Lo_D.size();
    }
    
    /* constructor passing wavelet name (only works for spln4, bior4.4, rbio2.8, and haar) */
    wavelet(int wavelet_name = SPLN4) {
        switch(wavelet_name) {
            case SPLN4:
                Lo_D = SPLN4_LO_D;
                Hi_D = SPLN4_HI_D;
                Lo_R = SPLN4_LO_R;
                Hi_R = SPLN4_HI_R;
                break;
            case BIOR4_4:
                Lo_D = BIOR4_4_LO_D;
                Hi_D = BIOR4_4_HI_D;
                Lo_R = BIOR4_4_LO_R;
                Hi_R = BIOR4_4_HI_R;
                break;
            case RBIO2_8:
                Lo_D = RBIO2_8_LO_D;
                Hi_D = RBIO2_8_HI_D;
                Lo_R = RBIO2_8_LO_R;
                Hi_R = RBIO2_8_HI_R;
                break;
            case HAAR:
                Lo_D = HAAR_LO_D;
                Hi_D = HAAR_HI_D;
                Lo_R = HAAR_LO_R;
                Hi_R = HAAR_HI_R;
                break;
            default:
                throw std::invalid_argument("wavelet::wavelet : invalid wavelet type"); 
        }
        length_filter = Lo_D.size();
    }
    
    /* n-level wavelet decomposition, each subband is a vector stored in 
       the output vector as {d1, d2, ..., dn, an}, being 1 the coarsest and
       n the finest level, thus having total size n+1 */
    std::vector<std::vector<T> > wavedec(const std::vector<T> &x, size_t n) const {
        #if SAFE_CHECK
        if(x.size()<(1<<n)) {
            throw std::invalid_argument("wavelet::wavedec : invalid decomposition size");
        }
        #endif
        std::vector<std::vector<T> > c(n+1);
        std::vector<T> a = x;
        for(size_t k=0; k<n; k++) {
            std::vector<T> ak;
            //          dk
            dwt(a, ak, c[k]);
            // save approximation coefficients for next iteration
            a = ak;
        }
        // store final approximation coefficients
        c[n] = a;
        return c;
    }
    
    /* DWT transform of 1D vector x, {a, d} coefficients are output arguments */
    void dwt(const std::vector<T> &x, std::vector<T> &a, std::vector<T> &d) const {
        size_t lx = x.size();
        size_t lenEXT = length_filter/2;
        size_t first = 1;
        size_t last = (size_t)(std::ceil((real_t)lx/2))*2;
        a = convdown(x, Lo_D, lenEXT, first, last);
        d = convdown(x, Hi_D, lenEXT, first, last);
    }
    
    /* n-level 2D wavelet decomposition, each subband is a matrix stored in
       the output vector as {d1, v1, h1, d2, v2, h2, ... dn, vn, hn, an}, 
       being 1 the coarsest and n the finest level, thus having size 3*n+1*/
    std::vector<matrix<T> > wavedec2(const matrix<T> &x, size_t n) const {
        #if SAFE_CHECK
        if(x.getM()<(1<<n) || x.getN()<(1<<n)) {
            throw std::invalid_argument("wavelet::wavedec2 : invalid decomposition size");
        }
        #endif
        std::vector<matrix<T> > c(3*n+1);
        matrix<T> a = x;
        for(size_t k=0; k<n; k++) {
            matrix<T> ak;
            //             hk        vk       dk
            dwt2(a, ak, c[3*k+2], c[3*k+1], c[3*k]);
            // save approximation coefficients for next iteration
            a = ak;
        }
        // store final approximation
        c[c.size()-1] = a;
        return c;
    }
    
    /* 2D DWT transform of matrix x, {a, h, v, d} coefficients are output arguments */
    void dwt2(const matrix<T> &x, matrix<T> &a, matrix<T> &h, matrix<T> &v, matrix<T> &d) const {
        size_t lenEXT = length_filter/2;
        std::vector<T> first = {1, 1};
        std::vector<T> last = {
            2*std::ceil((real_t)x.getM()/2), 
            2*std::ceil((real_t)x.getN()/2)
        };
        convdown2(LL_LH, x, lenEXT, first, last, a, h);
        convdown2(HL_HH, x, lenEXT, first, last, v, d);
    }
    
    /* calculate 2D subbands given matrix x:
       if mode=LL_LH => Lo_D applied to the rows => a and h subbands are obtained
       if mode=HL_HH => Hi_D applied to the rows => v and d subbands are obtained */
    void convdown2(int mode, const matrix<T> &x, size_t lenEXT, const std::vector<T> &first, 
            const std::vector<T> &last, matrix<T> &xL, matrix<T> &xH) const {
        
        // convolve and decimate each row with the chosen filter
        matrix<T> z = matrix<T>(x.getM(), (size_t)std::floor(last[1]/2));
        for(size_t k=0; k<x.getM(); k++) {
            std::vector<T> y = convdown(x.get_row(k), (mode==LL_LH? Lo_D : Hi_D), lenEXT, first[1], last[1]);
            z.set_row(y,k);
        }
        // obtain subband L/H H
        matrix<T> tmp = matrix<T>((size_t)std::floor(last[0]/2), (size_t)std::floor(last[1]/2));
        for(size_t k=0; k<z.getN(); k++) {
            std::vector<T> y = convdown(z.get_col(k), Hi_D, lenEXT, first[0], last[0]);
            tmp.set_col(y,k);
        }
        xH = tmp;
        // L/H L
        for(size_t k=0; k<z.getN(); k++) {
            std::vector<T> y = convdown(z.get_col(k), Lo_D, lenEXT, first[0], last[0]);
            tmp.set_col(y,k);
        }
        xL = tmp;
    }
    
    
    
    /* reconstruct signal of length sz from wavelet coeff c */
    std::vector<T> waverec(const std::vector<std::vector<T> > &c, size_t sz) const {
        return appcoef(c, 0, sz);
    }
    
    /* reconstruct 1D approximation of vector at level n, i.e. n=0 means full 
       reconstruction, n=1 means reconstruction at the second most 
       coarser level, and so on; c contains the d-level decomposition of the
       vecor (thus c has size L=d+1), then n<=d [=L-1]; the size of the
       reconstructed data is inferred directly by c unless n=0 (i.e., there
       is no greater level to infer the size from the coefficients vector c */
    std::vector<T> appcoef(const std::vector<std::vector<T> > &c, size_t n) const {
        #if SAFE_CHECK
        if(n==0 || n>(c.size()-1)) {
            throw std::invalid_argument("wavelet::appcoef : invalid reconstruction level");
        }
        #endif
        return appcoef(c, n, c.at(n-1).size());
    }
    
    /* reconstruct approximation of signal at level n having length sz */
    std::vector<T> appcoef(const std::vector<std::vector<T> > &c, size_t n, size_t sz) const {
        // c = {d1 d2 ... dn an}
        std::vector<T> a = c.back();
        for(size_t k=c.size()-1; k!=n; k--) {
            // get size from next upper level, if current level is the 
            // final then resort to the final lenght argument sz
            size_t s = k!=1? c.at(k-2).size() : sz;
            // reconstruct signal of size s from kth level approx and detail
            a = idwt( a, c[k-1], s );
        }
        return a;
    }
    
    /* inverse DWT transform {a(k),d(k)} -> a(k+1) of length s */
    std::vector<T> idwt(const std::vector<T> &a, const std::vector<T> &d, size_t s) const {
        return upconv(a,Lo_R,s) + upconv(d,Hi_R,s);
    }
    
    /* reconstruct 2D matrix of size sz[0]-by-sz[1] from wavelet coeff c */
    matrix<T> waverec2(const std::vector<matrix<T> > &c, const std::vector<size_t> &sz) const {
        // call appcoef with maximum level of reconstruction (0)
        return appcoef2(c, 0, sz);
    }
    
    /* reconstruct 2D approximation of matrix at level n, i.e. n=0 means full 
       reconstruction, n=1 means reconstruction at the second most 
       finer level; assuming c contains the d-level decomposition of a 2D
       signal (thus c has size L=3*d+1), then n<=d [=(L-1)/3]; the size of the
       reconstructed data is inferred directly by c; in case of n=0 there
       is no greater level to infer the size, thus the size of the
       reconstructed data is double the size of the second most finers level */
    matrix<T> appcoef2(const std::vector<matrix<T> > &c, size_t n) const {
        #if SAFE_CHECK
        if(n>(c.size()-1)/3) {
            throw std::invalid_argument("wavelet::appcoef2 : invalid reconstruction level");
        }
        #endif
        std::vector<size_t> sz = (n==0)? c[0].size()*(size_t)2 : c.at(3*n-1).size();
        return appcoef2(c, n, sz);
    }
    
    /* reconstruct 2D approximation of matrix at level n having size sz[0]-by-sz[1] */
    matrix<T> appcoef2(const std::vector<matrix<T> > &c, size_t n, 
            const std::vector<size_t> &sz) const {
        // c = {d1, v1, h1, d2, v2, h2, ... dn, vn, hn, an}
        matrix<T> a = c.back();
        for(size_t k=c.size()-1; k>3*n; k-=3) {
            // select size of reconstruction from any next upper subband, 
            // we use subband h (arbitrary as h, v, d have same size)
            std::vector<size_t> s = k!=3? c.at(k-4).size() : sz;
            //              hk      vk      dk
            a = idwt2( a, c[k-1], c[k-2], c[k-3], s );
        }
        return a;
    }
    
    /* get specific detail subband at level n from 2D coefficients c */
    matrix<T> detcoef2(const std::vector<matrix<T> > &c, int subband, size_t n) const {
        // c = {d1, v1, h1, d2, v2, h2, ... dn, vn, hn, an}
        #if SAFE_CHECK
        if(subband!=HORIZONTAL && subband!=VERTICAL && subband!=DIAGONAL) {
            throw std::invalid_argument("wavelet::detcoef2 : invalid subband");
        }
        if(n==0 || n>(c.size()-1)/3) {
            throw std::invalid_argument("wavelet::detcoef2 : invalid decomposition level");
        }
        #endif
        return c.at(3*(n-1)+subband);
    }
    
    static void set_appcoef2(std::vector<matrix<T> > &c, const matrix<T> &app) {
        #if SAFE_CHECK
        if(c.back().size()!=app.size()) {
            throw std::invalid_argument("wavelet::set_appcoef2 : invalid subband size");
        }
        #endif
        c.back() = app;
    }
    
    static void set_detcoef2(std::vector<matrix<T> > &c, const matrix<T> &band,
            int subband, size_t n) {
        // c = {d1, v1, h1, d2, v2, h2, ... dn, vn, hn, an}
        #if SAFE_CHECK
        if(subband!=HORIZONTAL && subband!=VERTICAL && subband!=DIAGONAL) {
            throw std::invalid_argument("wavelet::set_detcoef2 : invalid subband");
        }
        if(n==0 || n>(c.size()-1)/3) {
            throw std::invalid_argument("wavelet::set_detcoef2 : invalid decomposition level");
        }
        if(c.at(3*(n-1)+subband).size()!=band.size()) {
            throw std::invalid_argument("wavelet::set_detcoef2 : invalid subband size");
        }
        #endif
        c.at(3*(n-1)+subband) = band;
    }
    
    /* inverse 2D DWT transform {a(k),h(k),v(k),d(k)} -> a(k+1) of size s[0]-by-s[1] */
    matrix<T> idwt2(const matrix<T> &a, const matrix<T> &h, const matrix<T> &v, 
            const matrix<T> &d, const std::vector<size_t> &s) const {
        // sum of upsampled subbands
        return upconv2(a,Lo_R,Lo_R,s) + upconv2(h,Hi_R,Lo_R,s) +
                upconv2(v,Lo_R,Hi_R,s) + upconv2(d,Hi_R,Hi_R,s);
    }
    
    /* upsampling with zero-padding and convolution of matrix */
    matrix<T> upconv2(const matrix<T> &x, const std::vector<T> &f1, 
            const std::vector<T> &f2, const std::vector<size_t> &s) const {
        // create tmp matrix upsampling each column
        matrix<T> tmp = matrix<T>(s[0], x.getN());
        for(size_t k=0; k<x.getN(); k++) {
            // reconstruct column of size s[0], i.e. M of output
            tmp.set_col( upconv(x.get_col(k),f1,s[0]) , k);
        }
        // create final approximation matrix by upsampling rows of tmp, 
        // having size defined by argument s = {M, N}
        matrix<T> y = matrix<T>(s);
        for(size_t k=0; k<y.getM(); k++) {
            // reconstruct row of size s[1], i.e. N of output
            y.set_row( upconv(tmp.get_row(k),f2,s[1]) , k);
        }
        return y;
    }
    
    
    
    /* convolution and decimation  */
    std::vector<T> convdown(const std::vector<T> &x, const std::vector<T> &f, 
            size_t lenEXT, size_t first, size_t last) const {
        std::vector<T> z = wextend(x, lenEXT);
        z = conv1(z, f, VALID);
        std::vector<T> y(last/2);
        for(size_t k=first; k<=last; k+=2) {
            // first is 1, thus k is odd integer, thus k/2 === floor(k/2)
            y[k/2] = z[k];
        }
        return y;
    }
    
    /* upsampling with zero-padding dyadup and convolution  */
    std::vector<T> upconv(const std::vector<T> &x, const std::vector<T> &f, size_t s) const {
        size_t lf = f.size();
        std::vector<T> y = dyadup(x);
        y = wextend(y, lf/2);
        y = conv1(y, f, FULL);
        y.erase(y.begin(), y.begin()+lf-1);
        y.erase(y.begin()+s, y.end());
        return y;
    }
    
    /* double the size of x interleaving zeros in alternate odd positions
       thus obtaining {x1, 0, x2, 0, x3, ... xn, 0} */
    std::vector<T> dyadup(const std::vector<T> &x) const {
        std::vector<T> y;
        y.reserve(2*x.size());
        for(size_t k=0; k<x.size(); k++) {
            y.push_back( x[k] );
            y.push_back( 0 );
        }
        return y;
    }
    
    /* periodized extension of x */
    std::vector<T> wextend(const std::vector<T> &x, size_t lenEXT) const {
        std::vector<T> temp = x;
        // add extra sample if signal is odd
        if(x.size()%2==1) {
            temp.push_back(temp.back());
        }
        // handle cases when x is shorter than lenEXT
        size_t rep = lenEXT / temp.size(); // (size_t)std::floor(lenEXT/temp.size());
        size_t len = lenEXT % temp.size();
        std::vector<T> y;
        y.reserve(2*(len+rep*temp.size())+temp.size());
        // copy last len elements at the beginning
        y.insert(y.begin(), temp.end()-len, temp.end());
        for(size_t k=0; k<2*rep+1; k++) {
            y.insert(y.end(), temp.begin(), temp.end());
        }
        // copy first len elements at the end
        y.insert(y.end(), temp.begin(), temp.begin()+len);
        return y;
    }
    
    /* linearly reconstruct from low-pass coefficients of n-level decomposition 
       wavedec: thus if input a has length N, the output will be N*(2^n) 
       NOTE: default value for wavelet length reconstruction: length l(k) at 
       level k is double the length at the next finer level 2*l(k-1), but
       it works only for signal with power-of-2 lengths, otherwise waverec 
       routine will fail because of inconsistencies between the detail 
       subbands and the reconstructed approx subbands from previous level 
       (wrongly assumed of double length) */
    std::vector<T> linrec(const std::vector<T> &a, size_t n) const {
        // allocate vectors for wavelet coefficients
        // c = {d1 d2 ... dn an}
        std::vector<std::vector<real_t> > c;
        c.reserve(n+1);
        size_t N = a.size();
        for(size_t i=0; i<n; i++) {
            // set ith level detail coefficients to zero
            c.push_back( std::vector<real_t>(N*(1<<(n-i-1)), 0) );
        }
        // set approx coefficients in last position (low-pass subband)
        c.push_back( a );
        // call appcoef setting length of reconstruction to N*2^n
        return appcoef(c, 0, N*(1<<n));
    }
    
    /* same as before, but for 2D signals, linear reconstruction from 
       M-by-N matrix to M*2^n-by-N*2^n matrix */
    matrix<T> linrec2(const matrix<T> &a, size_t n) const {
        // allocate vector of matrices for wavelet coefficients
        // c = {d1 v1 h2 d2 v2 h2 ... dn vn hn an}
        std::vector<matrix<real_t> > c;
        c.reserve(3*n+1);
        size_t M = a.getM();
        size_t N = a.getN();
        for(size_t i=0; i<n; i++) {
            // set ith level diagonal, vertical, and horizontal coefficients to zero
            matrix<real_t> tmp = matrix<real_t>(M*(1<<(n-i-1)), N*(1<<(n-i-1)), 0);
            c.push_back( tmp );
            c.push_back( tmp );
            c.push_back( tmp );
        }
        // set approx coefficients in last position (low-pass subband)
        c.push_back( a );
        // call appcoef setting length of reconstruction to M*2^n-by-N*2^n
        return appcoef2(c, 0, {M*(1<<n),N*(1<<n)});
    }
    
    /* n-level linear upsampling of vector x from length N to length N*2^n
       using the reconstruction filter Lo_R */
    std::vector<T> linups(std::vector<T> x, size_t n) const {
        size_t lenEXT = length_filter/2; // !!! round(length_filter/2) but length_filter is always even
        /* double the size with zero-elems, apply padding, and convolve with Lo_R */
        for(size_t l=0; l<n; l++) {
            x = dyadup(x);
            x = wextend(x,lenEXT);
            x.pop_back();
            x = conv1(x, Lo_R, VALID);
        }
        return x;
    }
    
};

#endif

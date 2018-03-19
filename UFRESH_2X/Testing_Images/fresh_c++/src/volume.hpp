
#ifndef VOLUME_H
#define VOLUME_H

#include <iostream>
#include <stdexcept>
#include <vector>

#include "master.h"
#include "utils.hpp"
#include "matrix.hpp"

template<class T>
class volume {
private:
    size_t M; // number of rows
    size_t N; // number of columns
    size_t P; // number of planes 
    
    std::vector<T> data; // column-major vectorized content of the volume (a-la Matlab)
    
public:
    
    /* create a volume object with _M rows, _N columns, and _P planes */
    volume(size_t _M, size_t _N, size_t _P, T init = 0) {
        M = _M;
        N = _N;
        P = _P;
        data = std::vector<T>(M*N*P, init);
    }
    
    /* create a volume object from vector sz = {_M, _N, _P} */
    volume(const std::vector<size_t> &sz, T init = 0) {
        #if SAFE_CHECK
        if(sz.size()!=3) {
            throw std::invalid_argument("volume::volume : invalid size vector");
        }
        #endif
        M = sz[0];
        N = sz[1];
        P = sz[2];
        data = std::vector<T>(M*N*P, init);
    }
    
    /* create a volue object of size _M-by-_N-by-_P with vectorized _data */
    volume(size_t _M, size_t _N, size_t _P, const std::vector<T> &_data) {
        #if SAFE_CHECK
        if(_M*_N*_P>_data.size()) {
            throw std::invalid_argument("volume::volume : size and data do not agree");
        }
        #endif
        M = _M;
        N = _N;
        P = _P;
        data = std::vector<T>(_data.begin(),_data.begin()+_M*_N*_P);
    }
    
    /* default empty constructor */
    volume() { }
    
    /* get number of rows */
    size_t getM() const {
        return M;
    }
    
    /* get number of columns */
    size_t getN() const {
        return N;
    }
    
    /* get number of planes */
    size_t getP() const {
        return P;
    }
    
    /* get number of elements */
    size_t numel() const {
        return M*N*P;
    }
    
    /* return size of matrix as vector */
    std::vector<size_t> size() const {
        return { M, N, P };
    }
    
    /* get linearized index from subscript */
    size_t get_ind(size_t m, size_t n, size_t p) const {
        return p*(M*N) + n*M + m;
    }
    
    /* get subscript index along M from linearized index idx */
    size_t get_subM(size_t idx) const {
        return (idx%(M*N))%M;
    }
    
    /* get subscript index along N from linearized index idx */
    size_t get_subN(size_t idx) const {
        return (idx%(M*N))/M;
    }
    
    /* get subscript index along P from linearized index idx */
    size_t get_subP(size_t idx) const {
        return idx/(M*N);
    }
    
    /* get data */
    std::vector<T> get_data() const {
        return data;
    }
    
    /* set data */
    void set_data(const std::vector<T> &_data) {
        #if SAFE_CHECK
        if(data.size()!=_data.size()) {
            throw std::invalid_argument("volume::set_data : size and data do not agree");
        }
        #endif
        data = _data;
    }
    
    /* get max element in volume */
    T max() const {
        return *std::max_element(data.begin(), data.end());
    }
    
    /* get min element in volume */
    T min() const {
        return *std::min_element(data.begin(), data.end());
    }
    
    /* i-th plane occupies positions p*M*N+(0:M*N-1) */
    matrix<T> get_plane(size_t p) const {
        #if SAFE_CHECK
        if(p>=P) {
            throw std::out_of_range("volume::get_plane : invalid index");
        }
        #endif
        return matrix<T>(M,N,std::vector<T>(data.begin()+p*M*N, data.begin()+(p+1)*M*N));
    }
    
    /* i-th plane occupies positions p*M*N+(0:M*N-1) */
    void set_plane(const matrix<T> &plane, size_t p) {
        #if SAFE_CHECK
        if(p>=P || M!=plane.getM() || N!=plane.getN()) {
            throw std::out_of_range("volume::set_plane : invalid assignment");
        }
        #endif
        size_t p_idx = p*M*N;
        for(size_t i=0; i<plane.numel(); i++) {
            data[p_idx+i] = plane(i);
        }
    }
    void set_plane(const std::vector<T> &plane, size_t p) {
        #if SAFE_CHECK
        if(p>=P || M*N!=plane.size()) {
            throw std::out_of_range("volume::set_plane : invalid assignment");
        }
        #endif
        size_t p_idx = p*M*N;
        for(size_t i=0; i<plane.size(); i++) {
            data[p_idx+i] = plane[i];
        }
    }
    
    /* get n-th column of p-th plane */
    std::vector<T> get_col(size_t n, size_t p) const {
        #if SAFE_CHECK
        if(n>=N || p>=P) {
            throw std::out_of_range("volume::get_col : invalid indices");
        }
        #endif
        return std::vector<T>(data.begin()+p*M*N+n*M, data.begin()+p*M*N+(n+1)*M);
    }
    
    /* set n-th column of p-th plane */
    void set_col(const std::vector<T> &col, size_t n, size_t p) {
        #if SAFE_CHECK
        if(n>=N || p>=P || col.size()!=M) {
            throw std::out_of_range("volume::set_col : invalid indices");
        }
        #endif
        std::copy(col.begin(), col.end(), data.begin()+p*M*N+n*M);
    }
    
    /* get m-th row of p-th plane */
    std::vector<T> get_row(size_t m, size_t p) const {
        #if SAFE_CHECK
        if(m>=M || p>=P) {
            throw std::out_of_range("volume::get_row : invalid indices");
        }
        #endif
        std::vector<T> row;
        row.reserve(N);
        size_t p_idx = p*M*N;
        for(size_t k=0; k<N; k++) {
            row.push_back( data.at(p_idx+k*M+m) );
        }
        return row;
    }
    
    /* set m-th row of p-th plane */
    void set_row(const std::vector<T> &row, size_t m, size_t p) {
        #if SAFE_CHECK
        if(m>=M || p>=P || row.size()!=N) {
            throw std::out_of_range("volume::set_row : invalid indices");
        }
        #endif
        size_t p_idx = p*M*N;
        for(size_t k=0; k<N; k++) {
            data.at(p_idx+k*M+m) = row.at(k);
        }
    }
    
    /* get vector corresponding to the values in (m,n) position of all planes */
    std::vector<T> get_xsec(size_t m, size_t n) const {
        #if SAFE_CHECK
        if(m>=M || n>=N) {
            throw std::out_of_range("volume::get_xsec : invalid indices");
        }
        #endif
        std::vector<T> xsec;
        xsec.reserve(P);
        size_t mn_idx = n*M+m;
        for(size_t k=0; k<P; k++) {
            // value at position (m,n) of k-th plane
            xsec.push_back( data[ k*M*N+mn_idx ] );
        }
        return xsec;
    }
    
    /* set values in (m,n) position of all planes equal to xsec */
    void set_xsec(const std::vector<T> &xsec, size_t m, size_t n) {
        #if SAFE_CHECK
        if(m>M || n>N || xsec.size()!=P) {
            throw std::out_of_range("volume::set_xsec : invalid indices");
        }
        #endif
        size_t mn_idx = n*M+m;
        for(size_t k=0; k<P; k++) {
            data[ k*M*N+mn_idx ] = xsec[k];
        }
    }
    
    /* extract vectorized 3-D patch of size _Mx_Nx_P from data having the
       voxel at linearized position idx as top-left-front corner */
    std::vector<T> extract_patch(size_t _M, size_t _N, size_t _P, 
            size_t idx) const {
        #if SAFE_CHECK
        if(_N>N || _M>M || _P>P || idx>M*N*(P-_P)+M*(N-_N)+(M-_M)) {
            throw std::out_of_range("volume::extract_patch : invalid index");
        }
        #endif
        std::vector<real_t> patch_3d;
        patch_3d.reserve(_M*_N*_P);
        // for each plane spanned by the 3-D patch
        for(size_t k=0; k<_P; k++) {
            // extract _M-by-_N block from k-th plane
            for(size_t i=0; i<_N; i++) {
                patch_3d.insert(patch_3d.end(), 
                        data.begin()+idx+k*M*N+i*M, 
                        data.begin()+idx+k*M*N+i*M+_M );
            }
        }
        return patch_3d;
    }
    
    /* extract vectorized 3-D patch of size _M-by-_N-by-_P composed by _P
       blocks _M-by-_N whose top-left corners are stored in linearized form 
       in the trajectory vector; the trajectory vector has thus size _P; 
       note that the indices are given with respect to the size of the volume */
    std::vector<T> extract_patch(size_t _M, size_t _N, 
            const std::vector<size_t> &trajectory) const {
        size_t _P = trajectory.size();
        std::vector<real_t> patch_3d;
        patch_3d.reserve(_M*_N*_P);
        // for each element in trajectory
        for(size_t k=0; k<_P; k++) {
            // extract linearized index of k-th _M-by-_N block
            size_t p_idx = trajectory[k];
            // store block to patch_3d
            for(size_t i=0; i<_N; i++) {
                patch_3d.insert(patch_3d.end(), 
                        data.begin()+p_idx+i*M, 
                        data.begin()+p_idx+i*M+_M );
            }
        }
        return patch_3d;
    }
    
    /* extract vectorized 3-D patch of size _M-by-_N-by-_P from data having
       the voxel at position (m,n,p) as top-left-front corner */
    std::vector<T> extract_patch(size_t _M, size_t _N, size_t _P, 
            size_t m, size_t n, size_t p) const {
        #if SAFE_CHECK
        if(m+_M>M || n+_N>N || p+_P>P) {
            throw std::out_of_range("volume::extract_patch : invalid indices");
        }
        #endif
        // internal call to the above linearized method
        return extract_patch(_M, _N, _P, get_ind(m,n,p));
    }
    
    /* add vectorized 3-D patch of size _M-by-_N-by-_P to data having the
       voxel at linearized position idx as top-left-front corner */
    void add_patch(const std::vector<T> &patch_3d, 
            size_t _M, size_t _N, size_t _P, size_t idx) {
        #if SAFE_CHECK
        if(patch_3d.size()!=_M*_N*_P) {
            throw std::invalid_argument("volume::add_patch : size and patch do not agree");
        }
        #endif
        // for each plane spanned by the 3-D patch
        for(size_t k=0; k<_P; k++) {
            // offset of k-th block in data
            size_t dk_idx = idx + k*M*N;
            // offset of k-th block in patch_3d
            size_t pk_idx = k*_M*_N;
            for(size_t i=0; i<_N; i++) {
                // offset of i-th column of k-th block in data 
                size_t di_idx = dk_idx + i*M;
                // offset of i-th column of k-th block in patch_3d 
                size_t pi_idx = pk_idx + i*_M;
                // updating _M-by-_N block in data
                for(size_t j=0; j<_M; j++) {
                    data.at(di_idx+j) += patch_3d.at(pi_idx+j);
                }
            }
        }
    }
    
    /* add vectorized 3-D patch of size _M-by-_N-by-_T, being _T the size
       of the trajectory vector which contains the top-left indices of all
       blocks of size _M-by-_N composing patch_3d */
    void add_patch(const std::vector<T> &patch_3d, size_t _M, size_t _N, 
            const std::vector<size_t> &trajectory) {
        size_t _T = trajectory.size();
        #if SAFE_CHECK
        if(patch_3d.size()!=_M*_N*_T) {
            throw std::invalid_argument("volume::add_patch : size and patch do not agree");
        }
        #endif
        for(size_t k=0; k<_T; k++) {
            // offset of k-th _M-by-_N block in data
            size_t dk_idx = trajectory[k];
            // offset of k-th _M-by-_N block in patch_3d
            size_t pk_idx = k*_M*_N;
            for(size_t i=0; i<_N; i++) {
                // offset of i-th column of k-th block in data 
                size_t di_idx = dk_idx + i*M;
                // offset of i-th column of k-th block in patch_3d 
                size_t pi_idx = pk_idx + i*_M;
                // updating block in data
                for(size_t j=0; j<_M; j++) {
                    data.at(di_idx+j) += patch_3d.at(pi_idx+j);
                }
            }
        }
    }
    
    /* add vectorized 3-D patch of size _M-by-_N-by-_P to data having the
       voxel at position (m,n,p) as top-left-front corner */
    void add_patch(const std::vector<T> &patch_3d, 
            size_t _M, size_t _N, size_t _P, 
            size_t m, size_t n, size_t p) {
        #if SAFE_CHECK
        if(m+_M>M || n+_N>N || p+_P>P) {
            throw std::out_of_range("volume::add_patch : invalid indices");
        }
        #endif
        // internal call to the above linearized method
        add_patch(patch_3d, _M, _N, _P, get_ind(m,n,p));
    }
    
    /* volume addition */
    volume<T> operator+(volume<T> x) const {
        #if SAFE_CHECK
        if(x.getM()!=M || x.getN()!=N || x.getP()!=P) {
            throw std::invalid_argument("volume::+ : invalid dimensions");
        }
        #endif
        x.set_data( x.get_data() + data );
        return x;
    }
    
    /* scalar-volume multiplication <this>*alpha */
    volume<T> operator*(T alpha) const {
        return volume<T>(M, N, P, data * alpha);
    }
    
    /* scalar-volume addition <this>*alpha */
    volume<T> operator+(T alpha) const {
        return volume<T>(M, N, P, data + alpha);
    }
    
    /* equal operator for assignement operation */
    void operator=(const volume<T> &x) {
        M = x.getM();
        N = x.getN();
        P = x.getP();
        data = x.get_data();
    }
    
    /* parenthesis operator (.,.,.) for element indexing and assignment */
    T& operator()(size_t m, size_t n, size_t p) {
        #if SAFE_CHECK
        if(m>=M || n>=N || p>=P) {
            throw std::out_of_range("volume::() : invalid indices");
        }
        #endif
        return data[ get_ind(m,n,p) ];
    }
    
    /* parenthesis operator (.) for linearized indexing and assignment */
    T& operator()(size_t idx) {
        #if SAFE_CHECK
        if(idx>=data.size()) {
            throw std::out_of_range("volume::() : invalid index");
        }
        #endif
        return data[idx];
    }
    
    /* parenthesis operator (.) for linearized indexing */
    T operator()(size_t idx) const {
        #if SAFE_CHECK
        if(idx>=data.size()) {
            throw std::out_of_range("volume::() : invalid index");
        }
        #endif
        return data[idx];
    }
    
    /* print content to standard output */
    void print() const {
        for(size_t p=0; p<P; p++) {
            get_plane(p).print();
        }
    }
    
};

#endif

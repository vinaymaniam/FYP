
#ifndef MATRIX_H
#define MATRIX_H

#include <complex>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <vector>

#include "master.h"
#include "utils.hpp"

#define ROW_MAJOR 0
#define COL_MAJOR 1

template<class T>
        class matrix {
private:
    size_t M; // number of rows
    size_t N; // number of columns
    
    std::vector<T> data; // column-major vectorized content of the matrix
    
public:
    /* create a matrix object with M rows and N columns with values = init */
    matrix(size_t _M, size_t _N, T init = 0) {
        M = _M;
        N = _N;
        data = std::vector<T>(M*N, init);
    }
    
    /* create a matrix object of size M-by-N and column-major vectorized _data */
    matrix(size_t _M, size_t _N, const std::vector<T> &_data) {
        #if SAFE_CHECK
        if(_M*_N!=_data.size()) {
            throw std::invalid_argument("matrix::matrix : size and data do not agree");
        }
        #endif
        M = _M;
        N = _N;
        data = _data;
    }
    
    /* create a matrix object from size vector {M, N} */
    matrix(const std::vector<size_t> &sz, T init = 0) {
        #if SAFE_CHECK
        if(sz.size()!=2) {
            throw std::invalid_argument("matrix::matrix : invalid size vector");
        }
        #endif
        M = sz[0];
        N = sz[1];
        data = std::vector<T>(M*N, init);
    }
    
    /* create a matrix object of size M-by-N and column-major vectorized _data */
    matrix(const std::vector<size_t> &sz, const std::vector<T> &_data) {
        #if SAFE_CHECK
        if(sz.size()!=2 || prod1(sz)!=_data.size()) {
            throw std::invalid_argument("matrix::matrix : invalid matrix dimension and data");
        }
        #endif
        M = sz[0];
        N = sz[1];
        data = _data;
    }
    
    /* default empty constructor */
    matrix() { }
    
    /* get number of rows */
    size_t getM() const {
        return M;
    }
    
    /* get number of columns */
    size_t getN() const {
        return N;
    }
    
    /* get subscript index along M from linearized index idx */
    size_t get_subM(size_t idx) const {
        #if SAFE_CHECK
        if(idx>=M*N) {
            throw std::invalid_argument("matrix::get_subM : invalid index");
        }
        #endif
        return idx%M;
    }
    
    /* get subscript index along N from linearized index idx */
    size_t get_subN(size_t idx) const {
        #if SAFE_CHECK
        if(idx>=M*N) {
            throw std::invalid_argument("matrix::get_subN : invalid index");
        }
        #endif
        return idx/M;
    }
    
    /* get linearized index from subscript */
    size_t get_ind(size_t m, size_t n) const {
        #if SAFE_CHECK
        if(m>=M || n>=N) {
            throw std::invalid_argument("matrix::get_ind : invalid coordinates");
        }
        #endif
        return n*M + m;
    }
    
    /* get number of elements */
    size_t numel() const {
        return M*N;
    }
    
    /* return size of matrix as vector */
    std::vector<size_t> size() const {
        return { M, N };
    }
    
    /* ith column occupies positions i*M+(0:M-1) */
    void set_col(const std::vector<T> &col, size_t i) {
        #if SAFE_CHECK
        if(i>=N || M!=col.size()) {
            throw std::invalid_argument("matrix::set_col : invalid assignment");
        }
        #endif
        std::copy(col.begin(), col.end(), data.begin()+M*i);
    }

    /* ith row occupies positions (0:N-1)*M+i */
    void set_row(const std::vector<T> &row, size_t i) {
        #if SAFE_CHECK
        if(i>=M || N!=row.size()) {
            throw std::invalid_argument("matrix::set_row : invalid assignment");
        }
        #endif
        for(size_t k=0; k<N; k++) {
            data.at(k*M+i) = row.at(k);
        }
    }
    
    /* populate matrix from the vector x (default is col-major) */
    void set_data(const std::vector<T> &_data, int mode = COL_MAJOR) {
        #if SAFE_CHECK
        if(M*N!=_data.size()) {
            throw std::invalid_argument("matrix::set_data : dimensions do not agree");
        }
        #endif
        if(mode==ROW_MAJOR) {
            for(size_t i=0; i<M; i++) {
                for(size_t k=0; k<N; k++) {
                    data.at(k*M+i) = _data.at(i*N+k);
                }
            }
        } else {
            data = _data;
        }
    }
    
    /* get ith column */
    std::vector<T> get_col(size_t i) const {
        #if SAFE_CHECK
        if(i>=N) {
            throw std::out_of_range("matrix::get_col : invalid index");
        }
        #endif
        return std::vector<T>(data.begin()+M*i, data.begin()+M*(i+1));
    }
    
    /* get ith row */
    std::vector<T> get_row(size_t i) const {
        #if SAFE_CHECK
        if(i>=M) {
            throw std::out_of_range("matrix::get_row : invalid index");
        }
        #endif
        std::vector<T> row;
        row.reserve(N);
        for(size_t k=0; k<N; k++) {
            row.push_back( data.at(k*M+i) );
        }
        return row;
    }
    
    /* get vectorized internal data */
    std::vector<T> get_data(int mode = COL_MAJOR) const {
        std::vector<T> y;
        if(mode==ROW_MAJOR) {
            y.reserve(M*N);
            for(size_t i=0; i<M; i++) {
                for(size_t k=0; k<N; k++) {
                    y.push_back( data.at(k*M+i) );
                }
            }
        } else {
            y = data;
        }
        return y;
    }
    
    /* get max element in matrix */
    T max() const {
        return *std::max_element(data.begin(), data.end());
    }
    
    /* get min element in matrix */
    T min() const {
        return *std::min_element(data.begin(), data.end());
    }
    
    /* delete ith column */
    void delete_col(size_t i) {
        #if SAFE_CHECK
        if(i>=N) {
            throw std::out_of_range("matrix::delete_col : invalid index");
        }
        #endif
        data.erase(data.begin()+M*i, data.begin()+M*(i+1));
        N--;
    }
    
    /* delete ith column */
    void delete_row(size_t i) {
        #if SAFE_CHECK
        if(i>=M) {
            throw std::out_of_range("matrix::delete_row : invalid index");
        }
        #endif
        for(size_t k=0; k<N; k++) {
            data.erase(data.begin()+k*(M-1)+i);
        }
        M--;
    }
    
    /* reshape <this> matrix from M-by-N to _M-by-_N */
    void reshape(size_t _M, size_t _N) {
        #if SAFE_CHECK
        if(_M*_N!=M*N || _M*_N!=data.size()) {
            throw std::out_of_range("matrix::reshape : the numner of element must not change");
        }
        #endif
        M = _M;
        N = _N;
    }
    
    /* extract a patch of size _M-by-_N from data, having top-left 
       corner in linearized index idx, i.e. idx=n*M+m, being (m,n) the 2D
       subscript indices, analogous to data(m:m+_M-1,n:n+_N-1); the sub 
       matrix is returned in column-major vectorized form */
    std::vector<T> extract_patch(size_t _M, size_t _N, size_t idx) const {
        #if SAFE_CHECK
        if(_N>N || _M>M || idx>M*(N-_N)+(M-_M)) {
            throw std::out_of_range("matrix::extract_patch : invalid dimensions");
        }
        #endif
        std::vector<T> patch;
        patch.reserve(_M*_N);
        for(size_t i=0; i<_N; i++) {
            patch.insert(patch.end(), data.begin()+idx+i*M, data.begin()+idx+i*M+_M );
        }
        return patch;
    }

    /* same as before but the index of the sub matrix top-left pixel is 
       given in subscript 2D form (m,n) */
    std::vector<T> extract_patch(size_t _M, size_t _N, size_t m, size_t n) const {
        #if SAFE_CHECK
        if(m+_M>M || n+_N>N) {
            throw std::out_of_range("matrix::extract_patch : invalid dimensions");
        }
        #endif
        // internally we still call the memory-efficient function after
        // converting the subscripts indices into linearized form idx=n*M+m
        return extract_patch(_M, _N, n*M+m);
    }
    
    /* add patch x of size _M-by_N (given in column-major vectorized form)
       to <this> matrix data at linearized index idx, i.e. idx=n*M+m, being 
       (m,n) the 2D subscript indices, analogous to the following op
       data(m:m+_M-1,n:n+_M-1) += x */
    void add_patch(const std::vector<T> &x, size_t _M, size_t _N, size_t idx) {
        #if SAFE_CHECK
        if(_N>N || _M>M || idx>M*(N-_N)+(M-_M)) {
            throw std::out_of_range("matrix::add_patch : invalid dimensions");
        }
        #endif
        for(size_t i=0; i<_N; i++) {
            size_t di_idx = idx + i*M;
            size_t xi_idx = i*_M;
            for(size_t j=0; j<_M; j++) {
                data[di_idx+j] += x[xi_idx+j];
            }
        }
    }
    
    /* same as before but the index of the sub matrix top-left pixel is 
       given in subscript 2D form (m,n) */
    void add_patch(const std::vector<T> &x, size_t _M, size_t _N, size_t m, size_t n) {
        #if SAFE_CHECK
        if(m+_M>M || n+_N>N) {
            throw std::out_of_range("matrix::add_patch : invalid dimensions");
        }
        #endif
        // internally we still call the memory-efficient function after
        // converting the subscripts indices into linearized form idx=n*M+m
        add_patch(x, _M, _N, n*M+m);
    }
    
    /* extract sub matrix of size _M-by-_N and top-left pixel (m,n) from
       matrix x, then add it to <this> matrix data, i.e. the following op
       data(m:m+_M-1,n:n+_M-1) += x(m:m+_M-1,n:n+_M-1) */
    void add_sub_matrix(const matrix<T> &x, size_t _M, size_t _N, size_t m, size_t n) {
        #if SAFE_CHECK
        if(m+_M>M || n+_N>N || m+_M>x.getM() || n+_N>x.getN()) {
            throw std::out_of_range("matrix::add_sub_matrix : invalid dimensions");
        }
        #endif
        size_t d_idx = n*M+m;
        size_t x_idx = n*x.getM()+m;
        for(size_t i=0; i<_N; i++) {
            size_t di_idx = d_idx + i*M;
            size_t xi_idx = x_idx + i*x.getM();
            for(size_t j=0; j<_M; j++) {
                data[di_idx+j] += x(xi_idx+j);
            }
        }
    }
    
    /* get diagonal values in a vector */
    std::vector<T> diag() const {
        std::vector<T> d;
        size_t diag_length = MIN(M,N);
        d.reserve(diag_length);
        for(size_t i=0; i<diag_length; i++) {
            d.push_back( data.at(i*M+i) );
        }
        return d;
    }
    
    /* element-wise multiplication, as A.*B in Matlab */
    matrix<T> dot_star(matrix<T> x) {
        #if SAFE_CHECK
        if(x.getM()!=M || x.getN()!=N) {
            throw std::invalid_argument("matrix::dot_star : invalid dimensions");
        }
        #endif
        x.set_data( data * x.get_data());
        return x;
    }
    
    /* element-wise division, as A./B in Matlab */
    matrix<T> dot_slash(matrix<T> x) {
        #if SAFE_CHECK
        if(x.getM()!=M || x.getN()!=N) {
            throw std::invalid_argument("matrix::dot_slash : invalid dimensions");
        }
        #endif
        x.set_data( data / x.get_data());
        return x;
    }
    
    /* element-wise rounding elements in matrix */
    void round() {
        std::transform(data.begin(), data.end(), data.begin(), [](const T i) {
            return std::round(i);
        });
    }
    
    /* element-wise clipping, as MIN(MAX(x,m),M) */
    void clip(T m, T M) {
        std::transform(data.begin(), data.end(), data.begin(), [&m,&M](const T i) {
            return MIN( MAX(i,m), M );
        });
    }
    
    /* parenthesis operator (.,.) for element indexing and assignment */
    T& operator()(size_t i, size_t j) {
        #if SAFE_CHECK
        if(i>=M || j>=N) {
            throw std::out_of_range("matrix::() : invalid indices");
        }
        #endif
        return data.at(j*M+i);
    }
    
    /* parenthesis operator (.,.) for element indexing */
    T operator()(size_t i, size_t j) const {
        #if SAFE_CHECK
        if(i>=M || j>=N) {
            throw std::out_of_range("matrix::() : invalid indices");
        }
        #endif
        return data.at(j*M+i);
    }
    
    /* parenthesis operator (.) for linear element indexing and assignment */
    T& operator()(size_t i) {
        #if SAFE_CHECK
        if(i>=M*N) {
            throw std::out_of_range("matrix::() : invalid index");
        }
        #endif
        return data.at(i);
    }
    
    /* parenthesis operator (.) for linear element indexing */
    T operator()(size_t i) const {
        #if SAFE_CHECK
        if(i>=M*N) {
            throw std::out_of_range("matrix::() : invalid index");
        }
        #endif
        return data.at(i);
    }
    
    /* parenthesis operator ({.}) for vectorized element indexing */
    std::vector<T> operator()(const std::vector<size_t> &idx) const {
        std::vector<T> y;
        y.reserve(idx.size());
        std::transform(idx.begin(), idx.end(), std::back_inserter(y),
                [this](size_t i) {
            return data.at(i);
        });
        return y;
    }
    
    /* copy input matrix to <this> */
    void operator=(const matrix<T> &x) {
        M = x.getM();
        N = x.getN();
        data = x.get_data();
    }
    
    /* matrix addition */
    matrix<T> operator+(matrix<T> x) const {
        #if SAFE_CHECK
        if(x.getM()!=M || x.getN()!=N) {
            throw std::invalid_argument("matrix::+ : invalid dimensions");
        }
        #endif
        x.set_data( x.get_data() + data );
        return x;
    }
    
    /* right scalar matrix addition <this>+alpha */
    matrix<T> operator+(T alpha) const {
        return matrix<T>(M, N, data + alpha);
    }
    
    /* right scalar matrix subtraction <this>-alpha */
    matrix<T> operator-(T alpha) const {
        return matrix<T>(M, N, data - alpha);
    }
    
    /* matrix subtraction <this>-x */
    matrix<T> operator-(matrix<T> x) const {
        #if SAFE_CHECK
        if(x.getM()!=M || x.getN()!=N) {
            throw std::invalid_argument("matrix::- : invalid dimensions");
        }
        #endif
        x.set_data( data - x.get_data() );
        return x;
    }
    
    /* scalar-matrix multiplication <this>*alpha */
    matrix<T> operator*(T alpha) const {
        return matrix<T>(M, N, data * alpha);
    }
    
    /* scalar-matrix division <this>/alpha */
    matrix<T> operator/(T alpha) const {
        return matrix<T>(M, N, data / alpha);
    }
    
    /* element-wise power operator <this>.^alpha */
    matrix<T> operator^(T alpha) const {
        return matrix<T>(M, N, data ^ alpha);
    }
    
    /* matrix negation -<this> */
    matrix<T> operator-() const {
        return matrix<T>(M, N, -data);
    }
    
    /* transposition !<this> (through overloaded ! unary operator) */
    matrix<T> operator!() const {
        matrix<T> y = matrix<T>(N, M);
        y.set_data(data, ROW_MAJOR);
        return y;
    }
    
    void print() const {
        for(size_t k=0; k<M; k++) {
            std::vector<T> x = get_row(k);
            std::copy(x.begin(), x.end(), std::ostream_iterator<T>(std::cout, " "));
            std::cout << std::endl;
        }
    }
};

#endif

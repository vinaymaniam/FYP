
#ifndef UTILS_H
#define UTILS_H

#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "master.h"
    
/* convert vector to complex with zero imaginary parts */
template<class T>
std::vector<std::complex<T> > real2complex(const std::vector<T> &x) {
    std::vector<complex_t> y;
    y.reserve(x.size());
    for(size_t i=0; i<x.size(); i++) {
        y.push_back( { x[i], 0 } );
    }
    return y;
}

/* construct vector of uniformly spaced element from n1 to n2 separated 
   by a given step (which by default is one), i.e. (n1:step:n2) */
template<class T, class U, class S>
std::vector<T> range(U n1, S n2, real_t step = 1) {
    #if SAFE_CHECK
    if(n1>n2) {
        throw std::invalid_argument("range : invalid range definition");
    }
    #endif
    std::vector<T> v;
    size_t n = (size_t)std::floor(std::real(n2-n1)/step)+1;
    v.reserve(n);
    std::generate_n(std::back_inserter(v), n, [&n1,step]() {
        U nn = n1;
        n1 += step;
        return nn;
    });
    return v;
}

/* overload operator+ for element-wise vector addition x1+x2 */
template<class T>
std::vector<T> operator+(std::vector<T> x1, const std::vector<T> &x2) {
    #if SAFE_CHECK
    if(x1.size()!=x2.size()) {
        throw std::invalid_argument("std::vector operator+ : invalid dimensions");
    }
    #endif
    // this is actually implemented as x1 += x2, so the first parameter x1
    // is passed by copy and then overwritten with the result
    std::transform(x1.begin(), x1.end(), x2.begin(),
            x1.begin(), std::plus<T>());
    return x1;
}

/* overload operator+= for element-wise vector addition x1+x2 */
template<class T>
std::vector<T>& operator+=(std::vector<T> &x1, const std::vector<T> &x2) {
    #if SAFE_CHECK
    if(x1.size()!=x2.size()) {
        throw std::invalid_argument("std::vector operator+ : invalid dimensions");
    }
    #endif
    // the first parameter x1 is passed by reference, modified and returned
    std::transform(x1.begin(), x1.end(), x2.begin(),
            x1.begin(), std::plus<T>());
    return x1;
}

/* overload operator- for element-wise vector subtraction x1-x2 */
template<class T>
std::vector<T> operator-(std::vector<T> x1, const std::vector<T> &x2) {
    #if SAFE_CHECK
    if(x1.size()!=x2.size()) {
        throw std::invalid_argument("std::vector operator- : invalid dimensions");
    }
    #endif
    std::transform(x1.begin(), x1.end(), x2.begin(),
            x1.begin(), std::minus<T>());
    return x1;
}

/* overload operator-= for element-wise vector subtraction x1-x2 */
template<class T>
std::vector<T>& operator-=(std::vector<T> &x1, const std::vector<T> &x2) {
    #if SAFE_CHECK
    if(x1.size()!=x2.size()) {
        throw std::invalid_argument("std::vector operator- : invalid dimensions");
    }
    #endif
    std::transform(x1.begin(), x1.end(), x2.begin(),
            x1.begin(), std::minus<T>());
    return x1;
}

/* overload operator* for element-wise vector multiplication x1.*x2 */
template<class T>
std::vector<T> operator*(std::vector<T> x1, const std::vector<T> &x2) {
    #if SAFE_CHECK
    if(x1.size()!=x2.size()) {
        throw std::invalid_argument("std::vector operator* : invalid dimensions");
    }
    #endif
    std::transform(x1.begin(), x1.end(), x2.begin(),
            x1.begin(), std::multiplies<T>());
    return x1;
}

/* overload operator*= for element-wise vector multiplication x1.*x2 */
template<class T>
std::vector<T>& operator*=(std::vector<T> &x1, const std::vector<T> &x2) {
    #if SAFE_CHECK
    if(x1.size()!=x2.size()) {
        throw std::invalid_argument("std::vector operator* : invalid dimensions");
    }
    #endif
    std::transform(x1.begin(), x1.end(), x2.begin(),
            x1.begin(), std::multiplies<T>());
    return x1;
}

/* overload operator/ for element-wise vector division x1./x2 */
template<class T>
std::vector<T> operator/(std::vector<T> x1, const std::vector<T> &x2) {
    #if SAFE_CHECK
    if(x1.size()!=x2.size()) {
        throw std::invalid_argument("std::vector operator/ : invalid dimensions");
    }
    #endif
    std::transform(x1.begin(), x1.end(), x2.begin(),
            x1.begin(), std::divides<T>());
    return x1;
}

/* overload operator/= for element-wise vector division x1./x2 */
template<class T>
std::vector<T>& operator/=(std::vector<T> &x1, const std::vector<T> &x2) {
    #if SAFE_CHECK
    if(x1.size()!=x2.size()) {
        throw std::invalid_argument("std::vector operator/= : invalid dimensions");
    }
    #endif
    std::transform(x1.begin(), x1.end(), x2.begin(),
            x1.begin(), std::divides<T>());
    return x1;
}

/* overload operator- for vector negation -x */
template<class T>
std::vector<T> operator-(std::vector<T> x) {
    std::transform(x.begin(), x.end(), x.begin(), std::negate<T>());
    return x;
}

/* overload operator+ for left scalar vector addition alpha+x */
template<class T>
std::vector<T> operator+(T alpha, std::vector<T> x) {
    std::transform(x.begin(), x.end(), x.begin(), [&alpha](const T i) {
        return alpha + i;
    });
    return x;
}

/* overload operator+ for right scalar vector addition x+alpha */
template<class T>
std::vector<T> operator+(std::vector<T> x, T alpha) {
    return alpha + x;
}

/* overload operator+ for left scalar vector subtraction alpha-x */
template<class T>
std::vector<T> operator-(T alpha, std::vector<T> x) {
    std::transform(x.begin(), x.end(), x.begin(), [&alpha](const T i) {
        return alpha - i;
    });
    return x;
}

/* overload operator+ for right scalar vector subtraction x-alpha */
template<class T>
std::vector<T> operator-(std::vector<T> x, T alpha) {
    std::transform(x.begin(), x.end(), x.begin(), [&alpha](const T i) {
        return i - alpha;
    });
    return x;
}

/* overload operator* for left scalar-vector multiplication alpha*x */
template<class T>
std::vector<T> operator*(T alpha, std::vector<T> x) {
    std::transform(x.begin(), x.end(), x.begin(), [&alpha](const T i) {
        return alpha * i;
    });
    return x;
}

/* overload operator* for right scalar-vector multiplication x*alpha */
template<class T>
std::vector<T> operator*(std::vector<T> x, T alpha) {
    return alpha * x;
}

/* overload operator/ for left scalar-vector division alpha/x */
template<class T>
std::vector<T> operator/(T alpha, std::vector<T> x) {
    std::transform(x.begin(), x.end(), x.begin(), [&alpha](const T i) {
        return alpha / i;
    });
    return x;
}

/* overload operator/ for right scalar-vector division x/alpha */
template<class T>
std::vector<T> operator/(std::vector<T> x, T alpha) {
    std::transform(x.begin(), x.end(), x.begin(), [&alpha](const T i) {
        return i / alpha;
    });
    return x;
}

/* overload operator^ for vector power x^alpha */
template<class T>
std::vector<T> operator^(std::vector<T> x, T alpha) {
    std::transform(x.begin(), x.end(), x.begin(), [&alpha](const T i) {
        return std::pow(i, alpha);
    });
    return x;
}

/* sum of elements in vector sum(x) */
template<class T>
T sum1(const std::vector<T> &x) {
    T val = 0;
    val = std::accumulate(x.begin(), x.end(), val);
    return val;
}

/* prod of elements in vector prod(x) */
template<class T>
T prod1(const std::vector<T> &x) {
    T val = 1;
    val = std::accumulate(x.begin(), x.end(), val, std::multiplies<int>());
    return val;
}

/* vectorized abs operator */
template<class T>
std::vector<T> abs1(std::vector<T> x) {
    std::transform(x.begin(), x.end(), x.begin(), [](const T i) {
        return std::abs(i);
    });
    return x;
}

/* vectorized round operator */
template<class T>
std::vector<T> round1(std::vector<T> x) {
    std::transform(x.begin(), x.end(), x.begin(), [](const T i) {
        return std::round(i);
    });
    return x;
}

/* print vector to std out */
template<class T>
void print1(const std::vector<T> &x) {
    std::copy(x.begin(), x.end(), std::ostream_iterator<T>(std::cout, " "));
    std::cout << std::endl;
}

#endif

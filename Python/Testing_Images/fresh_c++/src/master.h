
#ifndef MASTER_H
#define MASTER_H

#include <chrono>
#include <complex>
#include <limits>

// define main numeric type
#define USE_DOUBLE 1 // boolean {0->float, 1->double}
// verbosity level
#define VERBOSE    1 // integer {0->no output, 1->some output, 2->more output, 3->lots of output}
// print warnings to standard output
#define WARNING    0 // boolean {0->don't show, 1->show}
// error checking
#define SAFE_CHECK 1 // boolean {0->loose error checking, 1->complete error checking}
// parallel computation
#ifndef PARALLEL
#define PARALLEL   0 // boolean {0->sequential, 1->parallel}
#endif

/******* do not modify under this line *******/

// define some common functions
#ifndef MAX
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif
#ifndef MIN
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#endif
#define FIX(N) ((N) < 0 ? std::ceil((N)) : std::floor((N)))
#define SQUARE(N) ((N)*(N))

#define DO_CAT(x,y) x ## y
#define CAT(x,y) DO_CAT(x,y)

// M_PI not defined under Win with VS C++ compiler
#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

// macros depending on main numeric type
#if USE_DOUBLE
typedef double real_t;
#define PREFIX_REAL d
#define PREFIX_CPLX z
#else
typedef float real_t;
#define PREFIX_REAL s
#define PREFIX_CPLX c
#endif
typedef std::complex<real_t> complex_t;

// smallest number that, when added to x, yields a result different from x,
// i.e., the machine epsilon (numeric precision) which depends on real_t
#define EPSILON std::numeric_limits<real_t>::epsilon() 
// define max value of size_t and real_t variables
#define MAX_SIZE_T std::numeric_limits<size_t>::max()
#define MAX_REAL_T std::numeric_limits<real_t>::max()

// chrono duration to milliseconds or micorseconds (for verbose mode)
#define MILLISEC(T) std::chrono::duration_cast<std::chrono::milliseconds>(T).count()
#define MICROSEC(T) std::chrono::duration_cast<std::chrono::microseconds>(T).count()

/* NOTES:
   Differences from Matlab implemetation:
    - upres is 8 for each iteration
    - indices of secondary (135deg) diagonal during FRI upsampling
    - pad is dependent on scale factor in shifted_resize
    - error correction is performed using exhaustive search block-matching
    - search windows in block-matching that exceeds image borders are
      cropped instead of being zero-padded
    - block matching is (optionally) performed before every linear mapping
    - number of discontinuities K is estimated from actual data
    - patch size  and eggs parameters in diagonal gradient fusion are fixed
    - linear mapping model M is computed via solving symmetric linear system
      from a Cholesky decomposition (without explicit matrix inversion)
    - shifted resize is (optionally) performed using cubic pchip interp1
      which allows to calculate extrapolated values
   Computational remarks:
    - MEX file compiled in MacOS (clang-700.1.76) returns exactly the same 
      results obtained by a pure Matlab implementation; MEX file compiled 
      under Win (VSC++2015) returns slightly different results; in both cases
      the same Intel MKL 11.1.1 BLAS/LAPACK libraries are used 
    - stand-alone binary compiled in MacOS using OpenBLAS 0.2.16 is unbearably 
      slow (issues in the compiler/OS/library combination?)
    - stand-alone binary compiled in MacOS using the Accelerate framework
      is roughly twice as slow as the MEX file using MKL (threading issues?)
      and suffers from numerical instabilities if float precision is used
 */

#endif

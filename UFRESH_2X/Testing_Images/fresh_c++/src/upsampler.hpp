
#ifndef UPSAMPLER_H
#define UPSAMPLER_H

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <vector>

#include "master.h"
#include "matrix.hpp"
#include "numeric.hpp"
#include "algebra.hpp"
#include "wavelet.hpp"
#include "fri.hpp"
#include "mapping.hpp"
#include "image.hpp"

// unmodifiable constants
#define STD_PROFILE    0
#define FAST_PROFILE   1
#define NO_DIAG_REG    0
#define DIAG_REG       1
#define NO_LIN_MAP     0
#define LIN_MAP        1
#define LO_RES_IN      0
#define HI_RES_IN      1
#define NO_FORCE_EVEN  0
#define FORCE_EVEN     1
#define DIAG_45        1
#define DIAG_135      -1

#define FRI_A          0
#define FRI_A_SYNTH    1
#define FRI_UPRES      2
#define FRI_SV_XY      3
#define FRI_SV_D       4
#define FRI_SHIFT      5
#define FRI_PSIZE      6
#define FRI_EGGS       7
#define LM_PSIZE       8
#define LM_SHIFT       9
#define LM_WIN         10
#define LM_K           11
#define LM_THRESHOLD   12
#define LM_FACTOR      13
#define LM_LAMBDA      14

// modifiable constants
#define MIN_LINE_LEN   0  // integer {0, 1, 2, ...} if 0 size of the image is used
#define DIAG_REG_ITER  1  // integer {1, 2, ...}
#define LIN_MAP_ITER   1  // integer {1, 2, ...}
#define REPEAT_BM      1  // boolean {0, 1}
#define SHI_RES_MODE   1  // integer {1->cubic, 2->pchip}
#define SHIFTED_RES    1  // boolean {0, 1}
#define MS_TIME_PAD    6  // integer {1, 2, ...}

class upsampler {
private:
    int filter;
    int filter_fri;
    real_t a;
    size_t upres;
    size_t log2_upres;
    
    int upsampling_mode;
    int correction_mode;
    int algorithm_mode;
    int profile;
    
    wavelet<real_t> wlt;
    wavelet<real_t> wlt_fri;
    std::unordered_map<size_t, fri> fri_map;
    
    std::unordered_map<int,real_t> params = {
        {FRI_A,        1.2},
        {FRI_A_SYNTH,  1},
        {FRI_UPRES,    8},
        {FRI_SV_XY,    0.3},
        {FRI_SV_D,     0.1},
        {FRI_SHIFT,    2},
        {FRI_PSIZE,    8},
        {FRI_EGGS,     10},
        {LM_PSIZE,     5},
        {LM_SHIFT,     1},
        {LM_WIN,       15},
        {LM_K,         4},
        {LM_THRESHOLD, 1000},
        {LM_FACTOR,    1.625},
        {LM_LAMBDA,    0.1}
    };
    
    #if PARALLEL==0 && VERBOSE>=3
    std::unordered_map<size_t, size_t> fri_counter_map;
    #endif
    
    #if PARALLEL
    static std::mutex mtx_fri_map;
    // number of concurrent threads supported by the architecture
    size_t cores = std::thread::hardware_concurrency();
    #endif
    
public:
    /* constructor */
    upsampler(int _filter, int _filter_fri, int _upsampling_mode = DIAG_REG, 
            int _correction_mode = LIN_MAP, int _algorithm_mode = LO_RES_IN,
            int _profile = STD_PROFILE) {
        #if SAFE_CHECK
        if(_upsampling_mode==DIAG_REG && DIAG_REG_ITER==0) {
            throw std::invalid_argument("upsampler::upsampler : invalid DIAG_REG_ITER");
        }
        if(_correction_mode==LIN_MAP && LIN_MAP_ITER==0) {
            throw std::invalid_argument("upsampler::upsampler : invalid LIN_MAP_ITER");
        }
        #endif
        // init internal parameters
        filter = _filter;
        filter_fri = _filter_fri;
        
        a = (algorithm_mode==LO_RES_IN)? get_param(FRI_A) : get_param(FRI_A_SYNTH);
        upres = (size_t)get_param(FRI_UPRES);
        log2_upres = (size_t)std::log2(upres);
        
        upsampling_mode = _upsampling_mode;
        correction_mode = _correction_mode;
        algorithm_mode = _algorithm_mode;
        profile = _profile;
        
        // init wavelet filters
        wlt = wavelet<real_t>(filter);
        wlt_fri = wavelet<real_t>(filter_fri);
        
        #if WARNING && MIN_LINE_LEN<16
        std::cout << "Warning : upsampler::upsampler : small MIN_LINE_LEN value" << std::endl;
        #endif
    }
    
    /* get parameter with id equal to key */
    real_t get_param(int key) const {
        if(params.find(key)==params.end()) {
            // if this point is reached, parameter key was not found
            throw std::invalid_argument("upsampler::get_param : invalid parameter key");
        }
        return params.at(key);
    }
    
    /* for each pM-by-pN block in I_ref, returns the locations of the most 
       similar K blocks within a search window in I_tar of size win-by-win 
       centered around the top-left pixel of the reference patch; similarity
       is computed as SSE (sum of squared-errors), and the reference patches
       are separated in each dimension by shift pixels; only patches having  
       SSE below or equal to threshold are considered similar; the function 
       returns a map structure with keys corresponding to the linearised  
       indices of the reference patch and values containing a std::vector of 
       size K containing the linearised indices of the best K patches (in 
       decresing order of similarity); if only k<K similar patches are found,
       the values in positions [k+1,K] in the vector are MAX_SIZE_T */
    std::unordered_map<size_t, std::vector<size_t> > block_matching(const matrix<real_t> &I_ref, 
            const matrix<real_t> &I_tar, size_t pM, size_t pN, size_t win, 
            size_t K, real_t threshhold, size_t shift) const {
        // compute ratio between I_ref and I_tar
        real_t ratioM = (real_t)I_ref.getM() / I_tar.getM();
        real_t ratioN = (real_t)I_ref.getN() / I_tar.getN();
        // win is centered around i, thus (i-half_win):(i+half_win)
        size_t half_win = win/2;
        std::unordered_map<size_t, std::vector<size_t> > best_K_idx_map;

        #if PARALLEL==0 && VERBOSE>=2
        std::chrono::high_resolution_clock::time_point t_start, t_end;
        t_start = std::chrono::high_resolution_clock::now();
        #endif
        
        for(size_t j=0; j<I_ref.getN()-pN+1; j++) {
            if(j%shift!=0 && j!=I_ref.getN()-pN) {
                continue;
            }
            // horizontal boundary conditions of search win
            size_t centerN = (size_t)std::round( (real_t)j/ratioN );
            size_t minN = (centerN<half_win)? 0 : centerN-half_win;
            size_t maxN = MIN(I_tar.getN()-pN+1, centerN+half_win+1);
            for(size_t i=0; i<I_ref.getM()-pM+1; i++) {
                if(i%shift!=0 && i!=I_ref.getM()-pM) {
                    continue;
                }
                // vertical boundary conditions of search win
                size_t centerM = (size_t)std::round( (real_t)i/ratioM );
                size_t minM = (centerM<half_win)? 0 : centerM-half_win;
                size_t maxM = MIN(I_tar.getM()-pM+1, centerM+half_win+1);
                // compute distance between patch_ref extracted from I_ref at
                // location (i,j) and every patch within a search window in
                // I_tar of size win-by-win centered around (centerM,centerN)
                std::vector<real_t> patch_ref = I_ref.extract_patch(pM, pN, i, j);
                // vector that contains the indices of the K patches within
                // the search window in I_tar centered at (centerM,centerN) 
                // having the K closest dst (sum of squared error) with  
                // respect to the reference patch extracted from I_ref
                std::vector<size_t> best_K_idx(K, MAX_SIZE_T);
                std::vector<real_t> best_K_dst(K, MAX_REAL_T);
                for(size_t n=minN; n<maxN; n++) {
                    for(size_t m=minM; m<maxM; m++) {
                        // extract target patch
                        std::vector<real_t> patch_tar = I_tar.extract_patch(pM, pN, m, n);
                        // compute patch distance as SSE (sum of squared errors)
                        real_t dst = squared_euclidean_distance(patch_ref, patch_tar);
                        // if current dst is smaller than the one relative 
                        // to the least similar patch AND than the threshold
                        if(dst<best_K_dst[K-1] && dst<threshhold) {
                            // update sse and idx in positions K-1
                            best_K_dst[K-1] = dst;
                            best_K_idx[K-1] = n*I_tar.getM()+m;
                            // find the position p in best_K_sse such that 
                            // best_K_dst[p-1]<=dst<best_K_dst[p], with 0<p<K
                            // and shift what is to the right by one, e.g.,
                            // with {s1, ..., si, sj, ..., sk} and p=j, we get
                            // {s1, ..., si, sp, sj, ...., sk-1}
                            for(size_t p=K-1; p!=0 && best_K_dst.at(p)<best_K_dst.at(p-1); p--) {
                                std::swap(best_K_dst[p],best_K_dst[p-1]);
                                std::swap(best_K_idx[p],best_K_idx[p-1]);
                            }
                        }
                    } // end for n
                } // end for m
                // add best K indices for patch (i,j) to output map
                best_K_idx_map[j*I_ref.getM()+i] = best_K_idx;
            } // end for i
        } // end for j

        #if PARALLEL==0 && VERBOSE>=2
        t_end = std::chrono::high_resolution_clock::now();
        std::cout << std::setw(MS_TIME_PAD) << MILLISEC(t_end-t_start) << 
                "ms      BL MATCH" << std::endl;
        #endif

        return best_K_idx_map;
    }

    /* applies linear mapping to correct each block pM-by-pN in I_rec
       using cross-scale similarities information contained in the K most 
       similar patches in I_l, whose indices are stored in best_K_idx_map; 
       I_l is the intermediate level image having size between the size of
       I_rec and its low-pass, whereas I_h is the "ground-truth" image
       corresponding to I_l, thus having the same size of I_l; the function
       extracts the K most similar patches from I_l and I_h, then learns the
       mapping between the patches in I_h (ground-truth) and I_l (low-res),
       then the same mapping is applied to correct the patches in I_rec */
    matrix<real_t> correction_by_linear_mapping(const matrix<real_t> &I_rec, 
            const matrix<real_t> &I_l, const matrix<real_t> &I_h, 
            const std::unordered_map<size_t,std::vector<size_t> > &best_K_idx_map,
            size_t pM, size_t pN) const {
        #if SAFE_CHECK
        if(I_l.getM()!=I_h.getM() || I_l.getN()!=I_h.getN()) {
            throw std::invalid_argument("upsampler::correction_by_linear_mapping : "
                    "dimensions of intermediate level images do not agree");
        }
        #endif
        // output variables
        matrix<real_t> buffer = matrix<real_t>(I_rec.size());
        matrix<real_t> weight = matrix<real_t>(I_rec.size());
        // internal parameters
        size_t s = pM*pN;
        real_t lambda = get_param(LM_LAMBDA);

        // weights used to aggergate error-corrected patches (uniform)
        std::vector<real_t> xw(s, 1);

        #if PARALLEL==0 && VERBOSE>=2
        std::chrono::high_resolution_clock::time_point t_start, t_end;
        t_start = std::chrono::high_resolution_clock::now();
        #endif

        // for each pair (idx -> {best_K_idx}) in map
        for(auto &pair : best_K_idx_map) {
            size_t idx = pair.first;
            std::vector<size_t> best_K_idx = pair.second;
            // calculate how many patches are found, i.e. n = sum(best_K_idx<Inf)
            size_t n = count_if(best_K_idx.begin(), best_K_idx.end(), [](const real_t i){
                return i<MAX_SIZE_T;
            });
            // low-resolution patch to correct
            std::vector<real_t> xl = I_rec.extract_patch(pM, pN, idx);
            // high-resolution corrected patch
            std::vector<real_t> xh;
            if(n==0) {
                // if no similar patches are found, then error correction
                // cannot be applied and thus xh=xl is not modified
                xh = xl;
            } else {
                // extract the n found patches from I_l and I_h and copy them
                // into the columns of Yl and Yh, respectively
                std::vector<real_t> Yl; Yl.reserve(s*n);
                std::vector<real_t> Yh; Yh.reserve(s*n);
                for(size_t k=0; k<n; k++) {
                    std::vector<real_t> xk = I_l.extract_patch(pM, pN, best_K_idx[k]);
                    Yl.insert(Yl.end(), xk.begin(), xk.end());
                    xk = I_h.extract_patch(pM, pN, best_K_idx[k]);
                    Yh.insert(Yh.end(), xk.begin(), xk.end());
                }
                // apply linear mapping on xl
                xh = linear_mapping(xl, Yl, Yh, s, n, lambda);
            }
            // update buffer and weight
            buffer.add_patch(xh, pM, pN, idx);
            weight.add_patch(xw, pM, pN, idx);
        }

        #if PARALLEL==0 && VERBOSE>=2
        t_end = std::chrono::high_resolution_clock::now();
        std::cout << std::setw(MS_TIME_PAD) << MILLISEC(t_end-t_start) << 
                "ms      ERR CORR" << std::endl;
        #endif

        // return weighted buffer
        return buffer.dot_slash(weight);
    }

    /* resize matrix m by a factor scale (as in function resize), then apply
       a shift proportional to factor to the resized result; the shift is 
       implemented via two cascading 1D cubic interp (first applied to the 
       rows and then to the columns) having 0 as extrapolation value; finally,
       the extrapolated zeros are then replaced with a reflective padding */
    matrix<real_t> shifted_resize(const matrix<real_t> &I, real_t scale,
            int shift_resize_mode = NO_FORCE_EVEN) const {
        return shifted_resize(I, scale, scale, shift_resize_mode);
    }
    matrix<real_t> shifted_resize(const matrix<real_t> &I,
            const std::vector<size_t> &dim, 
            int shift_resize_mode = NO_FORCE_EVEN) const {
        #if SAFE_CHECK
        if(dim.size()!=2) {
            throw std::invalid_argument("upsampler::shifted_resize : invalid dim vector");
        }
        #endif
        real_t scaleM = (real_t)dim[0] / I.getM();
        real_t scaleN = (real_t)dim[1] / I.getN();
        return shifted_resize(I, scaleM, scaleN, shift_resize_mode);
    }
    matrix<real_t> shifted_resize(const matrix<real_t> &I, real_t scaleM,
            real_t scaleN, int shift_resize_mode) const {
        #if SAFE_CHECK
        if(scaleM<=0 || scaleN<=0) {
            throw std::invalid_argument("upsampler::shifted_resize : invalid scale parameters");
        }
        #endif
        // resize matrix
        if(shift_resize_mode==FORCE_EVEN) {
            size_t Ms = (size_t)std::ceil(I.getM()*scaleM);
            size_t Ns = (size_t)std::ceil(I.getN()*scaleN);
            scaleM = (real_t)(Ms+(Ms%2==1))/I.getM();
            scaleN = (real_t)(Ns+(Ns%2==1))/I.getN();
        }
        matrix<real_t> I_res = resize(I, scaleM, scaleN, BICUBIC);
        #if SHIFTED_RES==0
        return I_res;
        #else
        size_t M = I_res.getM();
        size_t N = I_res.getN();
        // interp2 of resized matrix on shifted x-values (two cascading interp1)
        real_t shift = (real_t)(0.5*(scaleN-1));
        #if SHI_RES_MODE==1
        real_t extrapval = (real_t)0.0;
        #endif
        // original sample locations
        std::vector<real_t> x = range<real_t>(1,(int)N);
        // sample locations to be interpolated
        std::vector<real_t> xi = x+shift;
        matrix<real_t> I_shift_res = matrix<real_t>(M,N);
        // interp1_cubic along each row
        for(size_t i=0; i<M; i++) {
            #if SHI_RES_MODE==1
            std::vector<real_t> yi = interp1_cubic(x,I_res.get_row(i),xi, extrapval);
            #elif SHI_RES_MODE==2
            std::vector<real_t> yi = interp1_pchip(x,I_res.get_row(i),xi);
            #else
            throw std::invalid_argument("upsampler::shifted_resize : invalid SHI_RES_MODE");
            #endif
            I_shift_res.set_row(yi, i);
        }
        // as before, but interp1_cubic is now applied along each column
        shift = (real_t)(0.5*(scaleM-1));
        size_t pad = (size_t)std::ceil(scaleM/2);
        x = range<real_t>(1,(int)M);
        xi = x+shift;
        for(size_t i=0; i<N-pad*(SHI_RES_MODE==1); i++) {
            #if SHI_RES_MODE==1
            std::vector<real_t> yi = interp1_cubic(x,I_shift_res.get_col(i),xi, extrapval);
            #elif SHI_RES_MODE==2
            std::vector<real_t> yi = interp1_pchip(x,I_shift_res.get_col(i),xi);
            #else
            throw std::invalid_argument("upsampler::shifted_resize : invalid SHI_RES_MODE");
            #endif
            I_shift_res.set_col(yi, i);
        }
        #if SHI_RES_MODE==1
        // reflective padding of size pad on rows and columns
        for(size_t i=0; i<pad; i++) {
            I_shift_res.set_row( I_shift_res.get_row(M-2*pad+i), M-1-i );
            I_shift_res.set_col( I_shift_res.get_col(N-2*pad+i), N-1-i );
        }
        #endif
        return I_shift_res;
        #endif
    }
    
    /* get FRI kernel of size N from internal map structure, and creates a 
       new one if a kernel of desired size N is not found */
    fri get_fri_kernel(size_t N) {
        #if PARALLEL
        std::lock_guard<std::mutex> lock(mtx_fri_map);
        #endif
        // FRI object
        fri FRI;
        // initialize a FRI kernel for desired length, if not already in map
        auto fri_map_iter = fri_map.find(N);
        // 1 level of decomposition (upsampling factor equal to 2)
        size_t level = 1;
        if(fri_map_iter==fri_map.end()) {
            // not found
            #if PARALLEL==0 && VERBOSE>=3
            std::chrono::high_resolution_clock::time_point t_start, t_end;
            t_start = std::chrono::high_resolution_clock::now();
            #endif
            fri_map.insert(std::make_pair( N, fri(N, level, a, filter_fri)) );
            #if PARALLEL==0 && VERBOSE>=3
            t_end = std::chrono::high_resolution_clock::now();
            std::cout << std::setw(MS_TIME_PAD) << MILLISEC(t_end-t_start) << 
                    "ms      KER(" << N << ")" << std::endl;
            fri_counter_map.insert(std::make_pair( N, 1 ) );
            #endif
            FRI = fri_map.at(N);
        } else {
            #if PARALLEL==0 && VERBOSE>=3
            fri_counter_map.at(N) += 1;
            #endif
            FRI = fri_map_iter->second;
        }
        return FRI;
    }
    
    /* fusion of diagonal coefficients based on gradient of the image */
    inline matrix<real_t> diagonal_fusion(const matrix<real_t> &I_fri, 
            const matrix<real_t> &coef_d45, const matrix<real_t> &coef_d135) {
        #if SAFE_CHECK
        if(coef_d45.getM()!=coef_d135.getM() || coef_d45.getN()!=coef_d135.getN() || 
                I_fri.getM()!=coef_d135.getM() || I_fri.getN()!=coef_d135.getN()) {
            throw std::invalid_argument("upsampler::diagonal_fusion : invalid image size");
        }
        #endif
        // diagonal 45/135 fusion based on gradient of I_fri
        matrix<real_t> Gmag = matrix<real_t>(I_fri.size());
        matrix<real_t> Gdir = matrix<real_t>(I_fri.size());
        gradient(I_fri, Gmag, Gdir);

        size_t shift = get_param(FRI_SHIFT);
        size_t psize = get_param(FRI_PSIZE);
        size_t eggs = get_param(FRI_EGGS);

        if(profile==FAST_PROFILE) {
            shift *= 2;
        }

        matrix<real_t> coef_d = matrix<real_t>(I_fri.size());
        matrix<real_t> weight = matrix<real_t>(I_fri.size());
        std::vector<real_t> patch_weight(psize*psize, 1);
        std::vector<real_t> ang = {-135, 45, -45, 135};

        for(size_t i=0; i<I_fri.getM()-psize+1; i++) {
            if(i%shift!=0 && i!=I_fri.getM()-psize) {
                continue;
            }
            for(size_t j=0; j<I_fri.getN()-psize+1; j++) {
                if(j%shift!=0 && j!=I_fri.getN()-psize) {
                    continue;
                }
                // extract patch psize-by-psize having (i,j) as top-left corner
                std::vector<real_t> mag = Gmag.extract_patch(psize, psize, i, j);
                std::vector<real_t> dir = Gdir.extract_patch(psize, psize, i, j);

                // sort indexes based on values of mag sorted in descending order
                std::vector<size_t> idx_mag = range<size_t>((size_t)0,mag.size()-1);
                std::sort(idx_mag.begin(), idx_mag.end(), [&mag](size_t i, size_t j) {
                    return mag[i]>mag[j];
                });

                // find the most common direction in current patch
                std::vector<real_t> id_vec(ang.size());
                std::vector<real_t> dang(ang.size());
                for(size_t e=0; e<eggs; e++) {
                    dang = abs1( dir[idx_mag[e]]-ang );
                    size_t idx = std::distance( dang.begin(), 
                            std::min_element(dang.begin(), dang.end()) );
                    id_vec[idx] += 1;
                }
                size_t idx = std::distance( id_vec.begin(), 
                        std::max_element(id_vec.begin(), id_vec.end()) );

                // if idx is {0,1} then main dir is {-135,45}, i.e. coef_d45
                // if idx is {2,3} then main dir is {135,-45}, i.e. coef_d135
                coef_d.add_sub_matrix( (idx<2)? coef_d45 : coef_d135, 
                        psize, psize, i, j);
                weight.add_patch(patch_weight, psize, psize, i, j);
            }
        }
        coef_d = coef_d.dot_slash(weight);
        return coef_d;
    }
    
    /* FRI diagonal regularization */
    inline matrix<real_t> regularize_diagonals(const matrix<real_t> &I_fri, 
            int mode = DIAG_45) {
        real_t sv_thresh = get_param(FRI_SV_D);
        // get indices of all diagonals in I_fri as a single 1-D vector idx
        // which is roughly a zig-zag scan (with direction depending on mode)
        std::vector<size_t> idx;
        idx.reserve(I_fri.numel());
        for(size_t i=0; i<I_fri.getM()+mode; i++) {
            for(size_t j=i; j<=I_fri.numel()-1; j+=I_fri.getM()+mode) {
                idx.push_back(j);
            }
        }
        // getting vectorized data from I_fri sorted following idx
        std::vector<real_t> data = I_fri(idx);
        // getting minimum length of chunks to process in data
        size_t min_line_len = MIN_LINE_LEN<=0? I_fri.getM() : MIN_LINE_LEN*2;
        
        std::vector<real_t> data_reg;
        data_reg.reserve(data.size());
        // getting number of chunks
        size_t n_chunks = (data.size()<min_line_len)? 1 : data.size()/min_line_len;
        for(size_t i=0; i<n_chunks; i++) {
            // extract i-th chunk from data of length >=min_line_len
            std::vector<real_t> x(
                    data.begin() + i*min_line_len,
                    data.begin() + ((i==n_chunks-1)? data.size() : (i+1)*min_line_len)
            );
            // apply 1-level wavelet decomposition
            std::vector<std::vector<real_t> > c = wlt_fri.wavedec(x, 1);
            // extract low-pass
            std::vector<real_t> x_low = c.back();
            // log2_upres-level linear wavelet reconstruction 
            std::vector<real_t> x_high = wlt_fri.linrec(x, log2_upres);
            // getting FRI kernel and perform FRI upsampling
            fri FRI = get_fri_kernel(x_low.size());
            std::vector<real_t> d_pos = FRI.discontinuities_location(x_low, 
                    x_high.size(), upres, sv_thresh);
            x_high = FRI.upsample_line(x_high, d_pos);
            // extract log2_upres-level approximation and store it to data_reg
            c = wlt_fri.wavedec(x_high, log2_upres);
            data_reg.insert(data_reg.end(), c.back().begin(), c.back().end());
        }
        
        // unroll zig-zag data_reg to output matrix coef_d
        matrix<real_t> coef_d = matrix<real_t>(I_fri.size());
        for(size_t i=0; i<idx.size(); i++) {
            coef_d(idx.at(i)) = data_reg.at(i);
        }
        
        return coef_d;
    }
    
    /* FRI upsampling */
    inline matrix<real_t> upsample_lines(const matrix<real_t> &I_low, 
            size_t MM, size_t NN, size_t M, size_t N, int mode = ROW_MAJOR) {
        // linear upsampling
        std::vector<real_t> data(M*N);
        for(size_t i=0; i<N; i++) {
            // 1-level linear upsampling of i-th row or column, the result
            // will have double the size of the input
            std::vector<real_t> line_up = 
                    wlt.linups( mode==ROW_MAJOR? I_low.get_col(i) : I_low.get_row(i), 1 );
            // storing upsampled line to data in an interleaved fashion using
            // stride equal to N, so that the vector data contains the 
            // equivalent of the linearized transposed matrix that has the
            // upsampled lines as secondary dimension
            for(size_t j=0; j<M; j++) {
                data[j*N+i] = line_up[j];
            }
        }
        // FRI upsampling along scanlines
        matrix<real_t> coef_xy = matrix<real_t>(MM, NN);
        size_t min_line_len = MIN_LINE_LEN<=0? N : MIN_LINE_LEN;
        real_t sv_thresh = get_param(FRI_SV_XY);
        
        std::vector<real_t> data_high;
        data_high.reserve(data.size()*2);
        // total number of chunks having length in [min_line_len, 2*min_line_len)
        size_t n_chunks = (data.size()<min_line_len)? 1 : data.size()/min_line_len;
        for(size_t i=0; i<n_chunks; i++) {
            // extracting i-th chunk
            std::vector<real_t> x_low(
                    data.begin() + i*min_line_len,
                    data.begin() + ((i==n_chunks-1)? data.size() : (i+1)*min_line_len)
            );
            // linear reconstruction of i-th chunk
            std::vector<real_t> x_high = wlt_fri.linrec(x_low, log2_upres+1);
            // getting FRI kernel
            fri FRI = get_fri_kernel(x_low.size());
            // FRI upsampling
            std::vector<real_t> d_pos = FRI.discontinuities_location(x_low, 
                    x_high.size(), upres, sv_thresh);
            x_high = FRI.upsample_line(x_high, d_pos);
            // wavelet decomposition of FRI upsampled result
            std::vector<std::vector<real_t> > c = wlt_fri.wavedec(x_high, log2_upres);
            // the approximation coefficients of the FRI upsampling is the result
            data_high.insert(data_high.end(), c.back().begin(), c.back().end());
        }
        // storing data_high in row or column-major format
        coef_xy.set_data( data_high, mode );
        
        return coef_xy;
    }

    /* FRI upsampling of a low-pass image I_low M-by-N to an image scaled 
       by a factor equal to 2, i.e. size(I_low) -> 2*size(I_low); 
       I_app is the low-resolution approximation */
    matrix<real_t> upsamplexy(const matrix<real_t> &I_app, const matrix<real_t> &I_low, 
            real_t min_val, real_t max_val) {
        // image dimension before and after upsaling
        size_t M  = I_low.getM();
        size_t N  = I_low.getN();
        size_t NN = 2*N;
        size_t MM = 2*M;
        
        // level for wavelet reconstruction
        size_t level_rec = std::log2( I_low.getM() / I_app.getM() ) + 1;

        #if PARALLEL==0 && VERBOSE>=2
        std::chrono::high_resolution_clock::time_point t_start, t_end;
        t_start = std::chrono::high_resolution_clock::now();
        #endif
        
        // FRI upsample rows
        matrix<real_t> coef_x = upsample_lines(I_low, MM, NN, MM, N, ROW_MAJOR);
        coef_x.clip(min_val, max_val);

        #if PARALLEL==0 && VERBOSE>=2
        t_end = std::chrono::high_resolution_clock::now();
        std::cout << std::setw(MS_TIME_PAD) << MILLISEC(t_end-t_start) << 
                "ms      UPS ROW" << std::endl;
        t_start = std::chrono::high_resolution_clock::now();
        #endif
        
        // FRI upsample columns
        matrix<real_t> coef_y = upsample_lines(I_low, MM, NN, NN, M, COL_MAJOR);
        coef_y.clip(min_val, max_val);

        #if PARALLEL==0 && VERBOSE>=2
        t_end = std::chrono::high_resolution_clock::now();
        std::cout << std::setw(MS_TIME_PAD) << MILLISEC(t_end-t_start) << 
                "ms      UPS COL" << std::endl;
        t_start = std::chrono::high_resolution_clock::now();
        #endif

        // level_rec decomposition of FRI-upsampled coefficients
        std::vector<matrix<real_t> > cx = wlt.wavedec2(coef_x, level_rec);
        std::vector<matrix<real_t> > c  = wlt.wavedec2(coef_y, level_rec);
        // create reconstruction wavelet coefficients, they must be organized as
        // {d1, v1, h1, d2, v2, h2, ..., dn, vn, hn, an}
        for(size_t i=1; i<=level_rec; i++) {
            // extract VERTICAL detail subbands from the coef_x at level i
            // and replace coefficients in the reconstruction coefficients
            wavelet<real_t>::set_detcoef2(c, wlt.detcoef2(cx, VERTICAL, i), VERTICAL, i);
            // NOTE: we arbitrarily keep the diagonal coefficients of the 
            //       decomposition of coef_y (we could as well keep those
            //       of coef_x or take some combination of both)
        }
        // replace approximation (low-pass) subband with input
        wavelet<real_t>::set_appcoef2(c, I_app);
        // FRI reconstruction and clipping
        matrix<real_t> I_fri = wlt.waverec2(c, coef_x.size());
        I_fri.clip(min_val, max_val);

        #if PARALLEL==0 && VERBOSE>=2
        t_end = std::chrono::high_resolution_clock::now();
        std::cout << std::setw(MS_TIME_PAD) << MILLISEC(t_end-t_start) << 
                "ms      WLT REC" << std::endl;
        t_start = std::chrono::high_resolution_clock::now();
        #endif

        // the x/y FRI reconstructed I_fri is corrupted by jagging along 
        // the diagonals; these artifacts can be alleviated by applying FRI 
        // downsampling and upsampling (regularization) along each primary 
        // and secondary diagonal of I_fri, thus obtaining two FRI regularized 
        // images, which are merged based on the gradient information of I_fri;  
        // this process can be repeated DIAG_REG_ITER (with DIAG_REG_ITER>=0), 
        // by updating I_fri at the end of every iteration (thus I_fri at n=0 
        // is equal to that obtained by x/y FRI, whereas at any n>0, I_fri 
        // is the result of the diagonal regularization at iteration n-1)
        for(size_t d=0; upsampling_mode==DIAG_REG && d<DIAG_REG_ITER; d++) {
            // diagonal 45deg regularization
            matrix<real_t> coef_d45 = regularize_diagonals(I_fri, DIAG_45);
            coef_d45.clip(min_val, max_val);

            #if PARALLEL==0 && VERBOSE>=2
            t_end = std::chrono::high_resolution_clock::now();
            std::cout << std::setw(MS_TIME_PAD) << MILLISEC(t_end-t_start) << 
                    "ms      REG MAIN DIAG" << std::endl;
            t_start = std::chrono::high_resolution_clock::now();
            #endif

            // diagonal 135deg regularization
            matrix<real_t> coef_d135 = regularize_diagonals(I_fri, DIAG_135);
            coef_d135.clip(min_val, max_val);

            #if PARALLEL==0 && VERBOSE>=2
            t_end = std::chrono::high_resolution_clock::now();
            std::cout << std::setw(MS_TIME_PAD) << MILLISEC(t_end-t_start) << 
                    "ms      REG SEC DIAG" << std::endl;
            t_start = std::chrono::high_resolution_clock::now();
            #endif

            // diagonal 45/135 fusion based on gradient of I_fri
            matrix<real_t> coef_d = diagonal_fusion(I_fri, coef_d45, coef_d135);
            coef_d.clip(min_val, max_val);

            #if PARALLEL==0 && VERBOSE>=2
            t_end = std::chrono::high_resolution_clock::now();
            std::cout << std::setw(MS_TIME_PAD) << MILLISEC(t_end-t_start) << 
                    "ms      DIAG REG" << std::endl;
            t_start = std::chrono::high_resolution_clock::now();
            #endif

            // wavelet decomposition of the fused coefficients
            c = wlt.wavedec2(coef_d, level_rec);
            // reconstruction replacing low-pass (approximation coefficients)
            wavelet<real_t>::set_appcoef2(c, I_app);
            I_fri = wlt.waverec2(c, {MM,NN});
            I_fri.clip(min_val, max_val);

            #if PARALLEL==0 && VERBOSE>=2
            t_end = std::chrono::high_resolution_clock::now();
            std::cout << std::setw(MS_TIME_PAD) << MILLISEC(t_end-t_start) << 
                    "ms      WLT REC" << std::endl;
            t_start = std::chrono::high_resolution_clock::now();
            #endif
        }

        return I_fri;
    }

    /* upsample image I by a factor of 2, thus the size(I) -> 2*size(I), where
       I_int is the original low-res image, and I is its currently upsampled 
       version (in the first iteration they are equivalent), i.e., we have
       2^level*size(I_int)=size(I); in practice this function should be 
       called in a loop, having level as loop counter, from 1 to the desired 
       upsampling: in order to get an upsampling factor of 2, level=1; to get 
       an upsampling factor of 4, level=2; and so on; in each iteration this 
       function should be called keeping I_int fixed and equal to the input
       image, and I the upsampled result obtained from the prev iteration */
    matrix<real_t> upsample_by_two(const matrix<real_t> &I_int, 
            const matrix<real_t> &I, size_t level) {
        // extract minimum and maximum value of I_int
        real_t min_val = I_int.min();
        real_t max_val = I_int.max();

        #if PARALLEL==0 && VERBOSE
        std::chrono::high_resolution_clock::time_point t_begin, t_start, t_end;
        t_start = std::chrono::high_resolution_clock::now();
        t_begin = t_start;
        #endif

        // FRI reconstruction of I: size(I) -> 2*size(I)
        // level-th wavelet approximation coefficients are the original 
        // low-pass input image I_int scaled by a factor proportional to
        // the levels of the wavelet decomposition
        matrix<real_t> I_app = I_int*(real_t)(1<<level);
        // I is the currently upsampled result, and it is used as input in
        // this iteration (level) scaled by 2 because it represents the 
        // 1-level wavelet approximation coefficients (i.e., I is assumed 
        // to be the low-pass of an unknown hi-res image having double its size)
        matrix<real_t> I_fri = upsamplexy(I_app, I*2, min_val, max_val);

        #if PARALLEL==0 && VERBOSE
        t_end = std::chrono::high_resolution_clock::now();
        std::cout << std::setw(MS_TIME_PAD) << MILLISEC(t_end-t_start) << "ms  FRI " << 
                I.getM() << "x" << I.getN() << std::endl;
        #endif

        // error correction by linear mapping using intermediate level images
        if(correction_mode==LIN_MAP) {
            // reading parameters from internal structure
            size_t p_size = (size_t)get_param(LM_PSIZE);
            size_t win = (size_t)get_param(LM_WIN);
            size_t K = (size_t)get_param(LM_K);
            real_t threshold = get_param(LM_THRESHOLD)*SQUARE(p_size);
            size_t shift = (size_t)get_param(LM_SHIFT);
            real_t factor = get_param(LM_FACTOR);
            if(profile==FAST_PROFILE) {
                shift *= 2;
                win *= 0.75;
            }

            #if PARALLEL==0 && VERBOSE
            t_start = std::chrono::high_resolution_clock::now();
            #endif
            
            // extract 1-level wavelet approx I_low of I of size(I)/2
            std::vector<matrix<real_t> > c = wlt_fri.wavedec2(I, 1);
            matrix<real_t> I_low = c.back();
            // FRI reconstruction of low-pass version of I: size(I)/2 -> size(I)
            matrix<real_t> I0 = upsamplexy(I_low, I_low, min_val, max_val);
            // if I.M is odd, then (because of the upsampling) I0.M will be 
            // 2*ceil(I.M/2) = I.M+1, thus we should delete one row
            if(I.getM()%2==1) {
                I0.delete_row(I0.getM()-1);
            }
            // likewise for odd values of I.N, we sould delete one column
            if(I.getN()%2==1) {
                I0.delete_col(I0.getN()-1);
            }

            #if PARALLEL==0 && VERBOSE
            t_end = std::chrono::high_resolution_clock::now();
            std::cout << std::setw(MS_TIME_PAD) << MILLISEC(t_end-t_start) << "ms  FRI " << 
                    I_low.getM() << "x" << I_low.getN() << std::endl;
            t_start = std::chrono::high_resolution_clock::now();
            #endif
            
            // create intermediate level images by upsampling I and the FRI 
            // reconstruction I0 of the low-pass version of I using 
            // bilinear interpolation with scale equal to factor
            
            matrix<real_t> I_inter = shifted_resize(I, factor); 
            
  
            //matrix<real_t> I_inter = imresize(Baby,0.8);
            matrix<real_t> I_inter_fri = shifted_resize(I0, I_inter.size(), NO_FORCE_EVEN);
            
            // block matching using the FRI upsampled version I_fri of the
            // input I as reference, and the intermediate level I_inter_fri
            // as target, the matching is done for each patch in I_fri
            std::unordered_map<size_t,std::vector<size_t> > best_K_idx_map = 
                    block_matching(I_fri, I_inter_fri, p_size, p_size, win, K, threshold, shift);
            // for each reference patch, extract the K most similar patches 
            // from I_inter_fri and I_inter, then calculate the linear mapping 
            // that transforms the patches in I_inter_fri to those in the 
            // "ground-truth" I_inter; the resulting image I_fri_lm contains
            // an average of the corrected patches in I_fri
            matrix<real_t> I_fri_lm = correction_by_linear_mapping(I_fri, 
                    I_inter_fri, I_inter, best_K_idx_map, p_size, p_size);
            I_fri_lm.clip(min_val, max_val);
            // in case of artificial downsampling we substitute the low-pass
            // with I_app (i.e., the input scaled by a factor proportional
            // to the level of the wavelet transform)
            if(algorithm_mode==HI_RES_IN) {
                // wavelet decomposition and reconstruction replacing low-pass,
                // i.e., the exact low-pass wavelet coefficients
                std::vector<matrix<real_t> > c = wlt.wavedec2(I_fri_lm, level);
                wavelet<real_t>::set_appcoef2(c, I_app );
                I_fri_lm = wlt.waverec2(c, I_fri_lm.size());
                I_fri_lm.clip(min_val, max_val);
            }

            #if PARALLEL==0 && VERBOSE
            t_end = std::chrono::high_resolution_clock::now();
            std::cout << std::setw(MS_TIME_PAD) << MILLISEC(t_end-t_start) << "ms  MAP " << 
                    I_fri.getM() << "x" << I_fri.getN() << std::endl;
            t_start = std::chrono::high_resolution_clock::now();
            #endif

            // the last correction by linear mapping can be iterated, each 
            // time updating the ground-truth intermediate level (i.e., the 
            // corrected image obtained in the previous iteration), and the 
            // FRI reconstruction of the same intermediate image; in this 
            // way a new estimate of the upscaled image is obtained in each 
            // iteration and is reused in the following one as starting point
            for(size_t l=0; l<LIN_MAP_ITER; l++) {
                // downscale error-corrected image so that we can use it in 
                // the next regularization step as new intermediate level, 
                // i.e., 2*size(I) -> size(I)*factor
                I_inter = shifted_resize(I_fri_lm, factor/2, FORCE_EVEN);
                // compute 1-level wavelet decomposition of I_inter
                c = wlt_fri.wavedec2(I_inter, 1);
                // extract low-pass approximation
                I_low = c.back();
                // FRI regularization of the low-pass of the intermediate 
                // level I_inter: size(I)*factor/2 -> size(I)*factor
                I_inter_fri = upsamplexy(I_low, I_low, min_val, max_val);

                #if PARALLEL==0 && VERBOSE
                t_end = std::chrono::high_resolution_clock::now();
                std::cout << std::setw(MS_TIME_PAD) << MILLISEC(t_end-t_start) << 
                        "ms  FRI " << I_low.getM() << "x" << I_low.getN() << std::endl;
                t_start = std::chrono::high_resolution_clock::now();
                #endif

                // error correction by linear mapping
                #if REPEAT_BM
                // as before, the reference image is I_fri but the target 
                // is the FRI-regularized version I_inter_fri of the new 
                // I_inter, having size size(I)*factor; computational time
                // can be saved by reusing the previous matching information 
                if(profile==STD_PROFILE) {
                    best_K_idx_map = block_matching(I_fri, I_inter_fri, 
                            p_size, p_size, win, K, threshold, shift);
                }
                #endif
                // linear mapping using I_inter_fri as intermediate level
                // and I_inter as ground-truth of the intermediate level to
                // update all blocks in I_fri
                I_fri_lm = correction_by_linear_mapping(I_fri, I_inter_fri, 
                        I_inter, best_K_idx_map, p_size, p_size);
                I_fri_lm.clip(min_val, max_val);
                // in case of artificial downsampling
                if(algorithm_mode==HI_RES_IN) {
                    // wavelet decomposition and reconstruction replacing 
                    // low-pass, i.e., the exact low-pass wavelet coefficients
                    std::vector<matrix<real_t> > c = wlt.wavedec2(I_fri_lm, level);
                    wavelet<real_t>::set_appcoef2(c, I_app);
                    I_fri_lm = wlt.waverec2(c, I_fri_lm.size());
                    I_fri_lm.clip(min_val, max_val);
                }

                #if PARALLEL==0 && VERBOSE
                t_end = std::chrono::high_resolution_clock::now();
                std::cout << std::setw(MS_TIME_PAD) << MILLISEC(t_end-t_start) << 
                        "ms  MAP " << I_fri.getM() << "x" << I_fri.getN() << std::endl;
                t_start = std::chrono::high_resolution_clock::now();
                #endif
            }
            I_fri = I_fri_lm;
            #if PARALLEL==0 && VERBOSE
            std::cout << "--------" << std::endl; 
            #endif
        }

        #if PARALLEL==0 && VERBOSE
        std::cout << std::setw(MS_TIME_PAD) << MILLISEC(t_end-t_begin) << 
                "ms  TOT " << I.getM() << "x" << I.getN() << "->" <<
                I_fri.getM() << "x" << I_fri.getN() << std::endl << std::endl;
        #endif
        
        #if PARALLEL==0 && VERBOSE>=3
        std::cout << std::setw(MS_TIME_PAD) << "LEN" << "   ";
        for(auto &j : fri_map.begin()->second.get_time_log()) {
            std::cout << std::setw(MS_TIME_PAD) << j.first << "    ";
        }
        std::cout << std::endl;
        for(auto &i : fri_map) {
            std::cout << std::setw(MS_TIME_PAD) << i.first << " ";
            for(auto &j : i.second.get_time_log()) {
                real_t avg_us = (real_t)j.second/fri_counter_map.at(i.first);
                std::cout << std::setw(MS_TIME_PAD) << 
                        (size_t)std::round(avg_us) << "us  ";
            }
            std::cout << std::setw(MS_TIME_PAD) << "x" + 
                    std::to_string(fri_counter_map.at(i.first)) << std::endl;
        }
        std::cout << std::endl;
        #endif

        return I_fri;
    }
    
    /* upsample input image I by a factor equal to 2^level */
    matrix<real_t> upsample(const matrix<real_t> &I, size_t level) {
        // starting point of the iterative upsampling is the luminance
        matrix<real_t> I_up = I;
        for(size_t n=1; n<=level; n++) {
            // at n-th level we obtain an image (2^n*M)-by-(2^n*N)
            I_up = upsample_by_two(I, I_up, n);
        }
        return I_up;
    }
    
    #if PARALLEL
    /* upsample input image I by a factor equal to 2^level using threads */
    matrix<real_t> upsample_parallel(const matrix<real_t> &I, size_t level) {
        // definition of block size ans shift between blocks
        size_t Nblk = 32;
        size_t shift = Nblk-Nblk*0.15;
        // if size of image is not big enough, then let process sequentially
        if(MAX(I.getM(),I.getN())<=Nblk) {
            return upsample(I, level);
        }
        // internal parameters
        real_t factor = (1<<level);
        std::mutex mtx_blk;
        
        // weights for block aggregation (2-D Gaussian window)
        std::vector<real_t> w = gausswin2d(factor*Nblk, 3.5);
        
        // output matrices
        matrix<real_t> buffer = matrix<real_t>(I.size()*(size_t)factor, 0);
        matrix<real_t> weight = matrix<real_t>(buffer.size(), 0);
        
        // definition of a thread poll of size equal to the number of cores
        std::vector<std::thread> t(cores);
        // definition of grid that splits evenly the image across threads
        // from the factorization of cores
        std::vector<size_t> f = factorize(cores);
        std::vector<size_t> grid = {1, 1};
        size_t min_delta = MAX_SIZE_T;
        for(size_t i=1; i<f.size()-1; i++) {
            size_t alpha = 1;
            for(size_t j=0; j<=i; j++) {
                alpha *= f.at(j);
            }
            size_t beta = 1;
            for(size_t j=i+1; j<f.size(); j++) {
                beta *= f.at(j);
            }
            size_t delta = MAX(alpha,beta)-MIN(alpha,beta);
            if(delta<min_delta) {
                grid = {alpha, beta};
                min_delta = delta;
            }
        }
        #if VERBOSE
        std::cout << "Spawned " << cores << 
                " threads processing blocks of size " <<
                Nblk << "x" << Nblk << " and overlap " <<
                (Nblk-shift) << std::endl;
        #endif
        // each thread processes a different portion of the image
        for(size_t id=0; id<t.size(); id++) {
            t[id] = std::thread([&w, &Nblk, &shift, &mtx_blk, &buffer, &weight,
                    this, &I, &level, &grid, &factor, id]() {
                size_t counter = 0;
                // definition of portion top-left and bottom-right 
                // coordinates based on the thread id
                size_t minY = (id%grid.at(0))*I.getM()/grid.at(0);
                size_t maxY = (id%grid.at(0)+1)*I.getM()/grid.at(0);
                size_t minX = (id/grid.at(0))*I.getN()/grid.at(0);
                size_t maxX = (id/grid.at(0)+1)*I.getN()/grid.at(0);
                // for each pixel in the portion (skipping those not
                // corresponding to shift and those generating a block
                // that would fall out of bounds)
                for(size_t i=minY; i<maxY; i++) {
                    if(i+Nblk>I.getM()) {
                        continue;
                    }
                    if((i-minY)%shift!=0 && i+Nblk!=I.getM()) {
                        continue;
                    }
                    for(size_t j=minX; j<maxX; j++) {
                        if(j+Nblk>I.getN()) {
                            continue;
                        }
                        if((j-minX)%shift!=0 && j+Nblk!=I.getN()) {
                            continue;
                        }
                        try {
                            // extract block
                            matrix<real_t> blk = matrix<real_t>(Nblk, Nblk,
                                    I.extract_patch(Nblk, Nblk, i, j) );
                            // upsample block
                            matrix<real_t> blk_up = upsample(blk, level);
                            
                            // aggregate weighted block after obtaining mutex
                            std::lock_guard<std::mutex> lock(mtx_blk);
                            buffer.add_patch(blk_up.get_data() * w, 
                                    blk_up.getM(), blk_up.getN(),
                                    factor*i, factor*j);
                            weight.add_patch(w, blk_up.getM(), blk_up.getN(),
                                    factor*i, factor*j);
                        } catch(std::exception &e) {
                            std::cout << e.what() << std::endl;
                        }
                        counter++;
                    }
                }
                #if VERBOSE
                std::lock_guard<std::mutex> lock(mtx_blk);
                std::cout << "Thread " << id << " processed " <<
                        counter << " blocks in image portion of size " <<
                        (maxY-minY) << "x" << (maxX-minX) << std::endl;
                #endif
            });
        }
        // wait for all threads to finish
        for(size_t id=0; id<t.size(); id++) {
            t[id].join();
        }
        // return aggregated data 
        return buffer.dot_slash( weight );
    }
    #endif

};

#if PARALLEL
std::mutex upsampler::mtx_fri_map;
#endif

#endif

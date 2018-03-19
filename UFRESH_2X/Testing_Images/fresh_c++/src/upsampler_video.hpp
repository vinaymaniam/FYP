
#ifndef UPSAMPLER_VIDEO_H
#define UPSAMPLER_VIDEO_H

#include <complex>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "master.h"
#include "algebra.hpp"
#include "fourier.hpp"
#include "mapping.hpp"
#include "motion.hpp"
#include "numeric.hpp"
#include "volume.hpp"
#include "wavelet.hpp"

#define KAISER_10 {0.0579822579513560,0.0832530063727552,0.1055531669261077,0.1221716514612482,0.1310399172885329, \
                   0.1310399172885329,0.1221716514612482,0.1055531669261077,0.0832530063727552,0.0579822579513560}
#define KAISER_9 {0.0648344247775422,0.0964569234753488,0.1233740684158834,0.1414367824270763,0.1477956018082986, \
                  0.1414367824270763,0.1233740684158834,0.0964569234753488,0.0648344247775422}
#define KAISER_8 {0.0735376887281223,0.1142140887928174,0.1469894384182885,0.1652587840607718, \
                  0.1652587840607718,0.1469894384182885,0.1142140887928174,0.0735376887281223}
#define KAISER_7 {0.0849690854635967,0.1391494662164391,0.1790343091320881,0.1936942783757522, \
                  0.1790343091320881,0.1391494662164391,0.0849690854635967}
#define KAISER_6 {0.1006747349028220,0.1761788369907856,0.2231464281063924, \
                  0.2231464281063924,0.1761788369907856,0.1006747349028220}
#define KAISER_5 {0.1236796411180537,0.2353512128364889,0.2819382920909148, \
                  0.2353512128364889,0.1236796411180537}
#define KAISER_4 {0.1609242290117612,0.3390757709882388,0.3390757709882388,0.1609242290117612}
#define KAISER_3 {0.2336675003192803,0.5326649993614395,0.2336675003192803}
#define KAISER_2 {0.5,0.5}
#define KAISER_1 {1.0}

#define KAISER { KAISER_1, KAISER_2, KAISER_3, KAISER_4, KAISER_5, KAISER_6, KAISER_7, KAISER_8, KAISER_9, KAISER_10 }

#define STEP 1

/* Create 3-D Kaiser window of size pM-by-pN-by-pT */
volume<real_t> create_3D_kaiser_window(size_t pM, size_t pN, size_t pT) {
    std::vector<std::vector<real_t> > kaiser_db = KAISER;
    volume<real_t> win(pM,pN,pT, 1);
    if(pM>kaiser_db.size() || pN>kaiser_db.size() || pT>kaiser_db.size()) {
        #if WARNING
        std::cout << "create_3D_kaiser_window : unsupported size : " << 
                "uniform weights will be used" << std::endl;
        #endif
        return win;
    }
    std::vector<real_t> kaiser_2D(pM*pN, 0);
    std::vector<real_t> kaiser_1D = kaiser_db.at(pT-1);
    gemm(REAL, 'N','N', pM, pN, 1, kaiser_db.at(pM-1), pM, kaiser_db.at(pN-1), 1, kaiser_2D, pM);
    for(size_t t=0; t<pT; t++) {
        win.set_plane(kaiser_2D*kaiser_1D.at(t), t);
    }
    return win;
}

/* Interpolate video along its motion field MF */
volume<real_t> interpolation_along_motion_field(const volume<size_t> &MF, 
        const volume<real_t> &video, size_t pM, size_t pN, size_t pT, 
        real_t factor, int mode = LINEAR, int force_even = 0) {
    #if SAFE_CHECK
    if(MF.getM()!=video.getM() || MF.getN()!=video.getN() || MF.getP()!=video.getP()) {
        throw std::invalid_argument("interpolation_along_motion_field : size of motion field and video do not agree");
    }
    #endif
    size_t stepM = STEP;
    size_t stepN = STEP;
    size_t stepT = 1;
    
    size_t intp_len = (size_t)std::ceil((real_t)video.getP()*factor);
    if(force_even && intp_len%2!=0) {
        intp_len--;
        factor = (real_t)intp_len/video.getP();
    }
    
    // initialize video_intp using standard constant-in-time interpolation
    real_t w_ini = 0.000001;
    volume<real_t> video_intp = volume<real_t>(video.getM(),video.getN(),intp_len);
    volume<real_t> weight = volume<real_t>(video_intp.size(), w_ini);
    std::vector<real_t> x = range<real_t>(0,video.getP()-1);
    std::vector<real_t> xi = range<real_t>(0,intp_len-1)/factor;
    for(size_t j=0; j<video.getN(); j++) {
        for(size_t i=0; i<video.getM(); i++) {
            std::vector<real_t> y = video.get_xsec(i,j);
            std::vector<real_t> yi = (mode==LINEAR)? interp1_linear(x,y,xi, 0.0) : interp1_pchip(x,y,xi);
            // padding by replication for values out of bounds
            for(size_t n=0; n<xi.size(); n++) {
                if(xi[n]>x.back()) {
                    yi[n] = y.back();
                }
            }
            video_intp.set_xsec(yi*w_ini, i,j);
        }
    }
    
    for(size_t k=0; k<video.getP(); k++) {
        if(k%stepT!=0 && k!=video.getP()-pT) {
            continue;
        }
        for(size_t j=0; j<(video.getN()-pN+1); j++) {
            if(j%stepN!=0 && j!=video.getN()-pN) {
                continue;
            }
            for(size_t i=0; i<(video.getM()-pM+1); i++) {
                if(i%stepM!=0 && i!=video.getM()-pM) {
                    continue;
                }
                // extract trajectory corresponding to reference coordinate
                std::vector<size_t> trajectory_ref = 
                        get_trajectory(MF, i, j, k, pT);
                size_t _T = trajectory_ref.size();
                if(_T==1) {
                    continue;
                }
                // extract motion compensated volume
                volume<real_t> Vr = volume<real_t>(pM,pN,_T, 
                        video.extract_patch(pM, pN, trajectory_ref) );
                // interpolate missing positions in reference trajectory
                // extract and interpolate spatial coordinates of trajectory
                std::vector<std::vector<real_t> > traj_intp = 
                        interpolate_trajectory(MF, i,j,k, pM, pN, _T, factor, LINEAR);
                std::vector<real_t> traj_m_intp = traj_intp.at(0);
                std::vector<real_t> traj_n_intp = traj_intp.at(1);
                std::vector<real_t> traj_t_intp = traj_intp.at(2);
                std::vector<real_t> xi = traj_intp.at(3);
                std::vector<real_t> x = traj_intp.at(4);
                
                // it can happen that there are no interpolation points
                // available for the given patch temporal span
                if(xi.size()==0) {
                    continue;
                }
                
                #if SAFE_CHECK
                if(!(traj_m_intp.size()==traj_n_intp.size() && traj_t_intp.size()==traj_n_intp.size() &&
                        traj_t_intp.size()==xi.size())) {
                    throw std::logic_error("interpolation_along_motion_field : interpolation grids do not agree");
                }
                if(x.size()!=Vr.getP()) {
                    throw std::logic_error("interpolation_along_motion_field : error in interpolation grid");
                }
                #endif
                
                // motion-based interpolation and aggregation
                volume<real_t> win = create_3D_kaiser_window(pM,pN,xi.size());
                for(size_t n=0; n<pN; n++) {
                    for(size_t m=0; m<pM; m++) {
                        std::vector<real_t> y = Vr.get_xsec(m,n);
                        std::vector<real_t> yi;
                        if(mode==LINEAR) {
                            yi = interp1_linear(x, y, xi, 0.0);
                        } else {
                            yi = interp1_pchip(x, y, xi);
                        }
                        for(size_t t=0; t<yi.size(); t++) {
                            video_intp(traj_m_intp.at(t)+m,traj_n_intp.at(t)+n,traj_t_intp.at(t)) += win(m,n,t)*yi[t];
                            weight(traj_m_intp.at(t)+m,traj_n_intp.at(t)+n,traj_t_intp.at(t)) += win(m,n,t);
                        }
                    }
                }
                
            }
        }
    }
    
    // convex combination (averaging)
    for(size_t i=0; i<weight.numel(); i++) {
        if(weight(i)>0) {
            video_intp(i) /= weight(i);
        }
    }
    
    return video_intp;
}

/* Get location of low-pass with factor f in Fourier spectrum of size dim */
inline size_t get_low_pass_bound(size_t dim, size_t factor) {
    return (size_t)std::floor( (real_t)(factor*dim)/2 ) - FIX((real_t)dim/2);
}

/* Get registration shift */
inline real_t get_shift(size_t dim, size_t translation) {
    if(translation+1>FIX((real_t)dim/2)) {
        return (real_t)translation-(real_t)dim;
    } else {
        return (real_t)translation;
    }
}

/* Upsample video by a factor 2 along its motion field MF */
volume<real_t> upsample_by_two_using_nl_registration(const volume<size_t> &MF, 
        const volume<real_t> &video_intp, const volume<real_t> &video, size_t K,
        size_t pM, size_t pN, size_t pT, size_t win, real_t threshold) {
    #if SAFE_CHECK
    if(video_intp.getP()!=video.getP()*2) {
        throw std::invalid_argument("upsample_by_two_using_nl_registration : invalid interpolated input");
    }
    if(MF.getM()!=video.getM() || MF.getN()!=video.getN() || MF.getP()!=video.getP()) {
        throw std::invalid_argument("upsample_by_two_using_nl_registration : size of motion field and video do not agree");
    }
    #endif
    // we compute the squared normalized difference, thus the threshold 
    // should be squared as well (avoids computing the expensive sqrt)
    real_t squared_threshold = SQUARE(threshold);
    // initialize internal parameters
    size_t stepM = STEP;
    size_t stepN = STEP;
    size_t stepT = 1;
    size_t half_win = win/2;
    
    fourier DFT;
    
    // upsampling in time by 2
    size_t factor = 2;
    
    // initialize video_reg as video_intp with very small weight
    real_t w_ini = 0.000001;
    volume<real_t> video_reg = video_intp*w_ini;
    volume<real_t> weight = volume<real_t>(video_reg.size(), w_ini);
    
    for(size_t k=0; k<(video.getP()); k++) {
        if(k%stepT!=0 && k!=video.getN()-pT) {
            continue;
        }
        // temporal boundary conditions of search win
        size_t minT = (k<half_win)? 0 : k-half_win;
        size_t maxT = MIN(video.getP(), k+half_win+1);
        for(size_t j=0; j<(video.getN()-pN+1); j++) {
            if(j%stepN!=0 && j!=video.getN()-pN) {
                continue;
            }
            // horizontal boundary conditions of search win
            size_t minN = (j<half_win)? 0 : j-half_win;
            size_t maxN = MIN(video.getN()-pN+1, j+half_win+1);
            // for voxel i in row j of frame k
            for(size_t i=0; i<(video.getM()-pM+1); i++) {
                if(i%stepM!=0 && i!=video.getM()-pM) {
                    continue;
                }
                // vertical boundary conditions of search win
                size_t minM = (i<half_win)? 0 : i-half_win;
                size_t maxM = MIN(video.getM()-pM+1, i+half_win+1);
                
                // extract motion compensated volume at (i,j,k)
                std::vector<size_t> trajectory_ref = 
                        get_trajectory(MF, i, j, k, pT);
                size_t _T = trajectory_ref.size();
                if(_T==1) {
                    continue;
                }
                std::vector<real_t> patch_ref = 
                        video.extract_patch(pM, pN, trajectory_ref);
                real_t threshold_volume = squared_threshold*patch_ref.size();
                
                // vectors containing the linearized indices of the K most
                // similar patches and their similarities, respectively
                std::vector<size_t> best_K_idx(K, MAX_SIZE_T);
                std::vector<real_t> best_K_dst(K, MAX_REAL_T);
                // most similar patch is the reference patch_ref (distance=0)
                best_K_idx[0] = video.get_ind(i,j,k);
                best_K_dst[0] = 0;
                // compute distance between patch_ref and all patches within
                // a search window win-by-win centered around 
                // (centerM,centerN,centerT) in video_inter_fri
                for(size_t t=minT; t<maxT; t++) {
                    for(size_t n=minN; n<maxN; n++) {
                        for(size_t m=minM; m<maxM; m++) {
                            if(i==m && j==n && t==k) {
                                continue;
                            }
                            // extract target patch
                            std::vector<real_t> patch_tar = 
                                    video.extract_patch(pM, pN, 
                                    get_trajectory(MF, m, n, t, pT));
                            if(patch_tar.size()<patch_ref.size()) {
                                continue;
                            }
                            // compute patch distance as SSE (sum of squared errors)
                            real_t dst = squared_euclidean_distance(patch_ref, patch_tar);
                            // if current dst is smaller than the one relative
                            // to the least similar patch AND than the threshold
                            if(dst<best_K_dst[K-1] && dst<threshold_volume) {
                                // update sse and idx in positions K-1
                                best_K_dst[K-1] = dst;
                                best_K_idx[K-1] = video.get_ind(m, n, t);
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
                        }
                    }
                }
                
                // number of matches
                size_t Kr = count_if(best_K_idx.begin(), best_K_idx.end(), [](const real_t i){
                    return i<MAX_SIZE_T;
                });
                
                // 3-D DFT of reference patch
                volume<complex_t> F_Vr = DFT.fftshift3( DFT.fft3( 
                        volume<complex_t>(pM,pN,_T, real2complex(patch_ref)) ) );
                std::vector<real_t> m_shift(Kr,0);
                std::vector<real_t> n_shift(Kr,0);
                std::vector<real_t> t_shift(Kr,0);
                std::vector<std::vector<real_t> > Vm(Kr);
                for(size_t p=0; p<Kr; p++) {
                    std::vector<size_t> trajectory = get_trajectory(MF, best_K_idx[p], _T);
                    Vm[p] = video.extract_patch( pM,pN,trajectory );
                    if(p==0) {
                        // first position is occupied by reference volume
                        m_shift[p] = 0;
                        n_shift[p] = 0;
                        t_shift[p] = 0;
                    } else {
                        // 3-D DFT of target patch
                        volume<complex_t> F_Vm = DFT.fftshift3( DFT.fft3( 
                                volume<complex_t>(pM,pN,_T, real2complex( Vm[p] )) ) );
                        
                        // cross-correlation in Fourier domain between Vr and Vm
                        volume<complex_t> CC = volume<complex_t>( 
                                F_Vr.getM(), 
                                F_Vr.getN(), 
                                F_Vr.getP()*factor );
                        
                        for(size_t n=0; n<F_Vm.numel(); n++) {
                            F_Vm(n) = F_Vr(n) * std::conj(F_Vm(n));
                        }
                        
                        CC.add_patch(F_Vm.get_data(), pM,pN,_T, 
                                get_low_pass_bound(pM,1), 
                                get_low_pass_bound(pN,1), 
                                get_low_pass_bound(_T,factor));
                        CC = DFT.ifftshift3( CC );
                        
                        // invert 3-D DFT of cross-correlation
                        CC = DFT.ifft3( CC );
                        
                        // get maximum value of cross-correlation
                        size_t idx_max_CC = 0;
                        real_t max_CC = CC(0).real();
                        for(size_t n=1; n<CC.numel(); n++) {
                            if(CC(n).real()>max_CC) {
                                idx_max_CC = n;
                                max_CC = CC(n).real();
                            }
                        }
                        
                        // map shift to high-resolution grid of upsampled video
                        m_shift[p] = get_shift( CC.getM(), CC.get_subM(idx_max_CC) );
                        n_shift[p] = get_shift( CC.getN(), CC.get_subN(idx_max_CC) );
                        t_shift[p] = get_shift( CC.getP(), CC.get_subP(idx_max_CC) );
                        
                    }
                }
                // extract and interpolate spatial coordinates of trajectory
                std::vector<std::vector<real_t> > traj_intp = 
                        interpolate_trajectory(MF, i,j,k, pM, pN, _T, factor, LINEAR);
                std::vector<real_t> traj_m_intp = traj_intp[0];
                std::vector<real_t> traj_n_intp = traj_intp[1];
                
                // aggregate each registered volume in high resolution grid 
                // with respect to the coordinates of the reference volume
                size_t minT = factor*video.get_subP(trajectory_ref[0]);
                size_t maxT = factor*video.get_subP(trajectory_ref.back());
                volume<real_t> win = create_3D_kaiser_window(pM,pN,_T);
                for(size_t p=0; p<Kr; p++) {
                    
                    // skip if there is no translation
                    if(p>0 && ((m_shift[p]==0 && n_shift[p]==0 && t_shift[p]==0) || t_shift[p]==0)) {
                        continue;
                    }
                    
                    volume<real_t> Vp = volume<real_t>(pM,pN,_T, Vm[p]);
                    
                    // aggregation weight
                    real_t a = std::sqrt(best_K_dst[p]/Vp.numel());
                    real_t w = std::exp( -0.5 * SQUARE(a) );
                    
                    // aggregation in higher-resolution grid, note that
                    // the registered volume could fall out of bounds
                    for(int t=0; t<_T; t++) {
                        // coordinate of the registered t-th block in the volume
                        int t_v = (int)factor*(k+t)+t_shift[p];
                        if(t_v<minT || t_v>maxT) {
                            continue;
                        }
                        for(int n=0; n<pN; n++) {
                            int n_v = (int)(traj_n_intp.at(factor*t+t_shift[p])+n)+n_shift[p];
                            if(n_v<0 || n_v>=video_reg.getN()) {
                                continue;
                            }
                            for(int m=0; m<pM; m++) {
                                int m_v = (int)(traj_m_intp.at(factor*t+t_shift[p])+m)+m_shift[p];
                                if(m_v<0 || m_v>=video_reg.getM()) {
                                    continue;
                                }
                                video_reg(m_v,n_v,t_v) += Vp(m,n,t)*w*win(m,n,t);
                                weight(m_v,n_v,t_v) += w*win(m,n,t);
                            }
                        }
                    }
                    
                }
            }
        }
    }
    
    // convex combination (averaging)
    for(size_t i=0; i<weight.numel(); i++) {
        if(weight(i)>0) {
            video_reg(i) = video_reg(i)/weight(i);
        }
    }
    
    return video_reg;
    
}

/* Apply linear mapping on video using motion-compensated patches */
volume<real_t> correction_by_linear_mapping_with_motion_field(const volume<size_t> &MF, 
        const volume<size_t> &MF_inter, const volume<real_t> &video, 
        const volume<real_t> &video_inter, const volume<real_t> &video_inter_est,
        size_t K, size_t pM, size_t pN, size_t pT, size_t win, real_t threshold) {
    
    #if SAFE_CHECK
    if(video_inter.getM()!=video_inter_est.getM() || video_inter.getN()!=video_inter_est.getN() ||
            video_inter.getP()!=video_inter_est.getP()) {
        throw std::invalid_argument
                ("correction_by_linear_mapping_with_motion_field : sizes of intermediate level videos do not agree");
    }
    if(MF.getM()!=video.getM() || MF.getN()!=video.getN() ||
            MF.getP()!=video.getP()) {
        throw std::invalid_argument
                ("correction_by_linear_mapping_with_motion_field : sizes of motion field and video do not agree");
    }
    if(MF_inter.getM()!=video_inter.getM() || MF_inter.getN()!=video_inter.getN() ||
            MF_inter.getP()!=video_inter.getP()) {
        throw std::invalid_argument
                ("correction_by_linear_mapping_with_motion_field : sizes of motion field and intermediate video do not agree");
    }
    #endif
    // we compute the squared normalized difference, thus the threshold 
    // should be squared as well (avoids computing the expensive sqrt)
    real_t squared_threshold = SQUARE(threshold);
    // initialize internal parameters
    real_t lambda = 0.1;
    size_t stepM = STEP;
    size_t stepN = STEP;
    size_t stepT = 1;
    real_t ratioM = (real_t)video.getM() / video_inter_est.getM();
    real_t ratioN = (real_t)video.getN() / video_inter_est.getN();
    real_t ratioT = (real_t)video.getP() / video_inter_est.getP();
    size_t half_win = win/2;
    volume<real_t> video_lm = volume<real_t>(video.size());
    volume<real_t> weight = volume<real_t>(video.size());
    
    for(size_t k=0; k<video.getP(); k++) {
        if(k%stepT!=0 && k!=video.getN()-pT) {
            continue;
        }
        // temporal boundary conditions of search win
        size_t centerT = (size_t)std::round( (real_t)k/ratioT );
        size_t minT = (centerT<half_win)? 0 : centerT-half_win;
        size_t maxT = MIN(video_inter_est.getP(), centerT+half_win+1);
        for(size_t j=0; j<video.getN()-pN+1; j++) {
            if(j%stepN!=0 && j!=video.getN()-pN) {
                continue;
            }
            // horizontal boundary conditions of search win
            size_t centerN = (size_t)std::round( (real_t)j/ratioN );
            size_t minN = (centerN<half_win)? 0 : centerN-half_win;
            size_t maxN = MIN(video_inter_est.getN()-pN+1, centerN+half_win+1);
            for(size_t i=0; i<video.getM()-pM+1; i++) {
                if(i%stepM!=0 && i!=video.getM()-pM) {
                    continue;
                }
                // vertical boundary conditions of search win
                size_t centerM = (size_t)std::round( (real_t)i/ratioM );
                size_t minM = (centerM<half_win)? 0 : centerM-half_win;
                size_t maxM = MIN(video_inter_est.getM()-pM+1, centerM+half_win+1);
                
                // extract reference volume from video at (i,j,k)
                std::vector<size_t> trajectory_ref = 
                        get_trajectory(MF, i, j, k, pT);
                size_t _T = trajectory_ref.size();
                std::vector<real_t> patch_ref = 
                        video.extract_patch(pM, pN, trajectory_ref);
                size_t s = patch_ref.size();
                #if SAFE_CHECK
                if(s<pM*pN) {
                    throw std::logic_error("invalid length of spatiotemporal volume");
                }
                #endif
                real_t threshold_volume = squared_threshold*s;
                // compute distance between patch_ref and all patches within
                // a search window win-by-win centered around 
                // (centerM,centerN,centerT) in video_inter_fri
                std::vector<size_t> best_K_idx(K, MAX_SIZE_T);
                std::vector<real_t> best_K_dst(K, MAX_REAL_T);
                for(size_t t=minT; t<maxT; t++) {
                    for(size_t n=minN; n<maxN; n++) {
                        for(size_t m=minM; m<maxM; m++) {
                            // extract target patch
                            std::vector<real_t> patch_tar = 
                                    video_inter_est.extract_patch(pM, pN, 
                                    get_trajectory(MF_inter, m, n, t, _T));
                            if(patch_tar.size()<patch_ref.size()) {
                                continue;
                            }
                            // compute patch distance as SSE (sum of squared errors)
                            real_t dst = squared_euclidean_distance(patch_ref, patch_tar);
                            // if current dst is smaller than the one relative
                            // to the least similar patch AND than the threshold
                            if(dst<best_K_dst[K-1] && dst<threshold_volume) {
                                // update sse and idx in positions K-1
                                best_K_dst[K-1] = dst;
                                best_K_idx[K-1] = video_inter_est.get_ind(m, n, t);
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
                        }
                    }
                }
                // linear mapping of patch_ref using best_K_idx in video_inter_fri
                // get number of matches
                size_t c = count_if(best_K_idx.begin(), best_K_idx.end(), [](const real_t i){
                    return i<MAX_SIZE_T;
                });
                std::vector<real_t> xh;
                volume<real_t> xw = create_3D_kaiser_window(pM,pN,_T);
                // if more than 0 matches
                if(c>0) {
                    std::vector<real_t> Yl;
                    Yl.reserve(s*c);
                    std::vector<real_t> Yh;
                    Yh.reserve(s*c);
                    // extract patches from "ground-truth" and estimated
                    // intermediate levels video_inter and video_inter_fri
                    for(size_t k=0; k<c; k++) {
                        std::vector<real_t> tmp = video_inter_est.extract_patch(pM, pN, 
                                get_trajectory(MF_inter,best_K_idx[k],_T));
                        Yl.insert(Yl.end(), tmp.begin(), tmp.begin()+s);
                        tmp = video_inter.extract_patch(pM, pN, 
                                get_trajectory(MF_inter,best_K_idx[k],_T));
                        Yh.insert(Yh.end(), tmp.begin(), tmp.begin()+s);
                    }
                    xh = linear_mapping(patch_ref, Yl, Yh, s, c, lambda);
                } else {
                    xh = patch_ref;
                }
                // update buffer and weight
                video_lm.add_patch(xh*xw.get_data(), pM, pN, trajectory_ref);
                weight.add_patch(xw.get_data(), pM, pN, trajectory_ref);
            }
        }
    }
    
    // convex combination (averaging)
    for(size_t i=0; i<weight.numel(); i++) {
        if(weight(i)>0) {
            video_lm(i) /= weight(i);
        }
    }
    
    return video_lm;
    
}
    
/* Apply linear mapping on video */
volume<real_t> correction_by_linear_mapping(const volume<real_t> &video, 
        const volume<real_t> &video_inter, const volume<real_t> &video_inter_est,
        size_t K, size_t pM, size_t pN, size_t pT, size_t win, real_t threshold) {
    #if SAFE_CHECK
    if(video_inter.getM()!=video_inter_est.getM() || video_inter.getN()!=video_inter_est.getN() ||
            video_inter.getP()!=video_inter_est.getP()) {
        throw std::invalid_argument
                ("correction_by_linear_mapping : sizes of intermediate level videos do not agree");
    }
    #endif
    // we compute the squared normalized difference, thus the threshold 
    // should be squared as well (avoids computing the expensive sqrt)
    real_t squared_threshold = SQUARE(threshold);
    // initialize internal parameters
    size_t stepM = STEP;
    size_t stepN = STEP;
    size_t stepT = 1;
    size_t s  = pM*pN*pT;
    real_t lambda = 0.1;
    size_t half_win = win/2;
    real_t ratioM = (real_t)video.getM() / video_inter_est.getM();
    real_t ratioN = (real_t)video.getN() / video_inter_est.getN();
    real_t ratioT = (real_t)video.getP() / video_inter_est.getP();
    real_t threshold_volume = squared_threshold*s;
    
    volume<real_t> video_lm = volume<real_t>(video.size(), 0);
    volume<real_t> weight = volume<real_t>(video.size(), 0);
    volume<real_t> xw = create_3D_kaiser_window(pM,pN,pT);
    // for each frame k
    for(size_t k=0; k<video.getP()-pT+1; k++) {
        if(k%stepT!=0 && k!=video.getN()-pT) {
            continue;
        }
        // temporal boundary conditions of search win
        size_t centerT = (size_t)std::round( (real_t)k/ratioT );
        size_t minT = (centerT<half_win)? 0 : centerT-half_win;
        size_t maxT = MIN(video_inter.getP()-pT+1, centerT+half_win+1);
        // for each row j of frame k
        for(size_t j=0; j<video.getN()-pN+1; j++) {
            if(j%stepN!=0 && j!=video.getN()-pN) {
                continue;
            }
            // horizontal boundary conditions of search win
            size_t centerN = (size_t)std::round( (real_t)j/ratioN );
            size_t minN = (centerN<half_win)? 0 : centerN-half_win;
            size_t maxN = MIN(video_inter_est.getN()-pN+1, centerN+half_win+1);
            // for voxel i in row j of frame k
            for(size_t i=0; i<video.getM()-pM+1; i++) {
                if(i%stepM!=0 && i!=video.getM()-pM) {
                    continue;
                }
                // vertical boundary conditions of search win
                size_t centerM = (size_t)std::round( (real_t)i/ratioM );
                size_t minM = (centerM<half_win)? 0 : centerM-half_win;
                size_t maxM = MIN(video_inter_est.getM()-pM+1, centerM+half_win+1);
                // reference coordinate is (i,j,k)
                std::vector<real_t> patch_ref = 
                        video.extract_patch(pM, pN, pT, i, j, k);
                // compute distance between patch_ref and all patches within
                // a search window win-by-win centered around 
                // (centerM,centerN,centerT) in video_inter_fri
                std::vector<size_t> best_K_idx(K, MAX_SIZE_T);
                std::vector<real_t> best_K_dst(K, MAX_REAL_T);
                for(size_t t=minT; t<maxT; t++) {
                    for(size_t n=minN; n<maxN; n++) {
                        for(size_t m=minM; m<maxM; m++) {
                            // extract target patch
                            std::vector<real_t> patch_tar = 
                                    video_inter_est.extract_patch(pM, pN, pT, m, n, t);
                            // compute patch distance as SSE (sum of squared errors)
                            real_t dst = squared_euclidean_distance(patch_ref, patch_tar);
                            // if current dst is smaller than the one relative
                            // to the least similar patch AND than the threshold
                            if(dst<best_K_dst[K-1] && dst<threshold_volume) {
                                // update sse and idx in positions K-1
                                best_K_dst[K-1] = dst;
                                best_K_idx[K-1] = video_inter_est.get_ind(m, n, t);
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
                        }
                    }
                }
                // linear mapping of patch_ref using best_K_idx in video_inter_fri
                // get number of matches
                size_t c = count_if(best_K_idx.begin(), best_K_idx.end(), [](const real_t i){
                    return i<MAX_SIZE_T;
                });
                std::vector<real_t> xh;
                // if more than 0 matches
                if(c>0) {
                    std::vector<real_t> Yl;
                    Yl.reserve(s*c);
                    std::vector<real_t> Yh;
                    Yh.reserve(s*c);
                    // extract patches from "ground-truth" and estimated
                    // intermediate levels video_inter and video_inter_est
                    for(size_t k=0; k<c; k++) {
                        std::vector<real_t> tmp = video_inter_est.extract_patch(pM, pN, pT, best_K_idx[k]);
                        Yl.insert(Yl.end(), tmp.begin(), tmp.end());
                        tmp = video_inter.extract_patch(pM, pN, pT, best_K_idx[k]);
                        Yh.insert(Yh.end(), tmp.begin(), tmp.end());
                    }
                    xh = linear_mapping(patch_ref, Yl, Yh, s, c, lambda);
                } else {
                    xh = patch_ref;
                }
                // update buffer and weight
                video_lm.add_patch(xh*xw.get_data(), pM, pN, pT, i, j, k);
                weight.add_patch(xw.get_data(), pM, pN, pT, i, j, k);
            }
        }
    }
    
    // convex combination (averaging)
    for(size_t i=0; i<weight.numel(); i++) {
        if(weight(i)>0) {
            video_lm(i) /= weight(i);
        }
    }
    
    return video_lm;
}

#endif

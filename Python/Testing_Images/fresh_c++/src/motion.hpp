
#ifndef MOTION_H
#define MOTION_H

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "master.h"
#include "numeric.hpp"
#include "volume.hpp"
#include "wavelet.hpp"

#define FULL_SEARCH         0
#define NO_MOTION           1

#define MULTI_SCALE         0
#define SINGLE_SCALE        1

#define MIN_LEVEL_SIZE      60
#define MIN_HALF_WIN_SIZE   4

#define FORWARD_TIME        1
#define INVERSE_TIME       -1

#define CENTRE_OF_HEXAGON   0
#define LARGE_HEXAGON_SIZE  7
#define SMALL_HEXAGON_SIZE  4
#define INIT_VALUE         -2

#define REGULARIZE_DST      1
#define GRADIENT_DST        1

#define LINEAR              0
#define CUBIC               1

/* the mask defines which point of the new hexagon should be calculated
   depending on the position of the previous best match, each row i
   corresponds to the point (large_hexagon_dy[i], large_hexagon_dx[i]) */
const std::vector<int> large_hexagon_mask = {1,  1,  1,  1,  1,  1,  1,  /* (0,0)   */
                                             0,  1,  1,  0,  0,  0,  1,  /* (-2,-1) */
                                             0,  1,  1,  1,  0,  0,  0,  /* (-2,1)  */
                                             0,  0,  1,  1,  1,  0,  0,  /* (0,2)   */
                                             0,  0,  0,  1,  1,  1,  0,  /* (2,1)   */
                                             0,  0,  0,  0,  1,  1,  1,  /* (2,-1)  */
                                             0,  1,  0,  0,  0,  1,  1}; /* (0,-2)  */
const std::vector<int> large_hexagon_dy   = {0, -2, -2,  0,  2,  2,  0};
const std::vector<int> large_hexagon_dx   = {0, -1,  1,  2,  1, -1, -2};
const std::vector<int> small_hexagon_dy   = {-1,  0,  1,  0};
const std::vector<int> small_hexagon_dx   = { 0,  1,  0, -1};

inline real_t coordinate_distance(const matrix<size_t> &motion_vectors, 
        size_t yR, size_t xR, size_t yT, size_t xT) {
    return ( SQUARE( (real_t)yR-(real_t)yT ) + SQUARE( (real_t)xR-(real_t)xT ) ) / 
            ( SQUARE(motion_vectors.getM())+SQUARE(motion_vectors.getN()) );
}

real_t get_median_direction(const matrix<size_t> &matching_table, 
        size_t win, size_t yR, size_t xR) {
    #if SAFE_CHECK
    if(yR>=matching_table.getM() || xR>=matching_table.getN()) {
        throw std::out_of_range("get_median_direction : invalid coordinates");
    }
    #endif
    // definition of the search window sround (yR,xR)
    size_t half_win = MAX(1,win/2);
    size_t minX = (xR<half_win)? 0 : xR-half_win;
    size_t maxX = MIN(matching_table.getN(), xR+half_win+1);
    size_t minY = (yR<half_win)? 0 : yR-half_win;
    size_t maxY = MIN(matching_table.getM(), yR+half_win+1);
    std::vector<real_t> dir_vec;
    dir_vec.reserve(win*win-1);
    for(size_t x=minX; x<maxX; x++) {
        for(size_t y=minY; y<maxY; y++) {
            // get coordinate of block that best matches the one at (y,x)
            size_t idx = matching_table(y,x);
            if(idx==MAX_SIZE_T) {
                continue;
            }
            real_t vy = (real_t)matching_table.get_subM(idx)-(real_t)y;
            real_t vx = (real_t)matching_table.get_subN(idx)-(real_t)x;
            dir_vec.push_back( (std::atan2(-vy,vx)+M_PI)/(2*M_PI) );
        }
    }
    real_t median_dir = -1;
    if(dir_vec.size()>0) {
        std::sort(dir_vec.begin(), dir_vec.end());
        size_t sz = dir_vec.size();
        median_dir = (sz%2==0)? (dir_vec[sz/2-1]+dir_vec[sz/2])/2.0 : dir_vec[sz/2];
    }
    return median_dir;
}

/* computes all motion vectors from reference frameR to target frameT and 
   returns a matrix of the same size of frameR containing in each position
   (i,j) the linearized index coordinate of the block in frameT which is
   the most similar to the reference block in frameR having its top-left
   pixel at position (i,j); the blocks have size pM-by-pN */
matrix<size_t> get_motion_vectors(const matrix<real_t> &frameR, const matrix<real_t> &frameT, 
        size_t offset, size_t pM, size_t pN, size_t half_win, real_t threshold, 
        int matching_mode = FULL_SEARCH, int estimation_mode = SINGLE_SCALE) {
    #if SAFE_CHECK
    if(frameR.getM()!=frameT.getM() || frameR.getN()!=frameT.getN()) {
        throw std::invalid_argument("get_motion_vectors : dimensions of frames do not agree");
    }
    #endif
    // we compute the squared normalized difference, thus the threshold 
    // should be squared as well (avoids computing the expensive sqrt)
    real_t squared_threshold = SQUARE(threshold);
    // window for computation of median direction
    size_t win_med_dir = 3;
    // definition of convex weights for distance metric
    real_t delta_1 = 0.85;
    real_t delta_2 = 0.10;
    real_t delta_3 = 0.05;
    // levels of multi-scale estimation
    size_t levels = 0;
    if(estimation_mode==MULTI_SCALE) {
        // get number of levels in multi-scale pyramid, thus we will have a
        // pyramid composed by frames having size scaled as 2^level, with 
        // level = {levels, levels-1, ... 0}, i.e.:
        // {frameR.size()/2^level, frameR.size()/2^(level-1), ..., frameR.size()/2^0}
        levels = (size_t)std::floor( std::log2((real_t)MIN(frameR.getM(),frameR.getN())/MIN_LEVEL_SIZE) );
    }
    // wavelet decomposition of input images
    wavelet<real_t> wlt(HAAR);
    std::vector<matrix<real_t> > coefR = wlt.wavedec2(frameR, levels);
    std::vector<matrix<real_t> > coefT = wlt.wavedec2(frameT, levels);
    // vector containing matrices of linearized indices; for each level in
    // the multiscale decomposition, we construct a matrix of the same size
    // of the image pairs, and then we populate such matrix putting in each
    // position (y,x) the linearized index of the block in target which
    // is the most similar to the reference block extracted at (y,x)
    std::vector<matrix<size_t> > matching_table_pyramid(levels+1);
    for(size_t level=0; level<=levels; level++) {
        // initializing useful parameters
        size_t L = levels-level;
        size_t scale = (1<<L);
        size_t half_win_L = MAX(MIN_HALF_WIN_SIZE,half_win/(1<<level));
        size_t pM_L = pM;
        size_t pN_L = pN;
        
        // offset to add to the linearized coordinate to be consistent
        // to the position of the input frames in the video
        size_t offset_L = (level==levels)? offset : 0;
        
        // extracting approximation at current level
        matrix<real_t> frameR_L = wlt.appcoef2(coefR, L) / (real_t)scale;
        matrix<real_t> frameT_L = wlt.appcoef2(coefT, L) / (real_t)scale;
        
        // initializing table at current level
        matching_table_pyramid[level] = matrix<size_t>(frameR_L.size(),MAX_SIZE_T);
        
        // for each pixel in the current level
        for(size_t xR=0; xR<frameR_L.getN()-pN_L+1; xR++) {
            for(size_t yR=0; yR<frameR_L.getM()-pM_L+1; yR++) {
                // position around which center the search window
                size_t yR_L = yR;
                size_t xR_L = xR;
                real_t median_dir = -1;
                // update the position if a lower level is available
                if(level>0) {
                    // coordinates corresponding to (yR,xR) in the lower level
                    size_t yR_low = yR/2;
                    size_t xR_low = xR/2;
                    // getting coordinate of the most similar block in the lower level
                    size_t idxT = matching_table_pyramid[level-1](yR_low,xR_low);
                    if(idxT!=MAX_SIZE_T) {
                        // if there is a match get coordinates of idx
                        size_t yT_low = matching_table_pyramid[level-1].get_subM(idxT);
                        size_t xT_low = matching_table_pyramid[level-1].get_subN(idxT);
                        // get the motion vector corresponding to the matched
                        // pair at the lower level; the motion vector
                        // should be multiplied by 2 in order to compensate
                        // for the change in scale from the lower level
                        int vy = 2*((int)yT_low-(int)yR_low);
                        int vx = 2*((int)xT_low-(int)xR_low);
                        // add the motion vector to the starting reference
                        // position at the current level
                        yR_L = yR + vy;
                        xR_L = xR + vx;
                        #if SAFE_CHECK
                        // check that the new coordinates do not fall out of bounds
                        if(yR_L>=frameR_L.numel() || xR_L>=frameR_L.numel()) {
                            throw std::out_of_range("get_motion_vectors : invalid reference coordinate");
                        }
                        #endif
                    }
                    if(delta_3>0) {
                        // compute median of direction within a neighbourhood
                        // centered around (yR_low,xR_low) in the lower scale
                        median_dir = get_median_direction(matching_table_pyramid[level-1],
                                win_med_dir, yR_low, xR_low);
                    }
                }
                // extract reference block from frameR_L
                std::vector<real_t> block_ref = frameR_L.extract_patch(pM_L, pN_L, yR, xR);
                // define a search window around the reference coordinate
                size_t minX = (xR_L<half_win_L)? 0 : xR_L-half_win_L;
                size_t maxX = MIN(frameT_L.getN()-pN_L+1, xR_L+half_win_L+1);
                size_t minY = (yR_L<half_win_L)? 0 : yR_L-half_win_L;
                size_t maxY = MIN(frameT_L.getM()-pM_L+1, yR_L+half_win_L+1);
                #if SAFE_CHECK
                if(minX>maxX || minY>maxY) {
                    throw std::logic_error("get_motion_vectors : invalid search window");
                }
                #endif
                // search for the block in frameT that best matches block_ref
                size_t best_block_idxT = MAX_SIZE_T;
                real_t min_l2norm = MAX_REAL_T;
                // for each position in search window
                for(size_t x=minX; x<maxX; x++) {
                    for(size_t y=minY; y<maxY; y++) {
                        // extract target block at position (y,x) in frameT_L
                        std::vector<real_t> block_tar = frameT_L.extract_patch(pM_L, pN_L, y, x);
                        // compute normalized l2 norm between ref block and target
                        real_t l2norm = delta_1 * squared_euclidean_distance(block_ref, block_tar) / (pM_L*pN_L);
                        // regularization term on the predicted position
                        if(delta_2>0) {
                            // incorporate delta_3 if there is no available
                            // median direction
                            real_t delta_23 = (median_dir>0)? delta_2 : delta_2 + delta_3;
                            real_t dy = (real_t)yR_L-(real_t)y;
                            real_t dx = (real_t)xR_L-(real_t)x;
                            l2norm += delta_23 * 0.5 * ( 
                                    SQUARE(dy)/(real_t)SQUARE(frameT_L.getM()) + 
                                    SQUARE(dx)/(real_t)SQUARE(frameT_L.getN()) );
                        }
                        // regularization term on the median direction
                        if(delta_3>0 && median_dir>0) {
                            real_t vy = (real_t)y-(real_t)yR;
                            real_t vx = (real_t)x-(real_t)xR;
                            real_t dir_diff = ( (std::atan2(-vy,vx)+M_PI)/(2*M_PI)-median_dir );
                            l2norm += delta_3 * SQUARE(dir_diff);
                        }
                        // if the distance is less than threshold and best so far
                        if(l2norm<min_l2norm && l2norm<squared_threshold) {
                            min_l2norm = l2norm;
                            best_block_idxT = frameR_L.get_ind(y, x);
                        }
                    }
                }
                // save matched position in table
                matching_table_pyramid[level](yR,xR) = best_block_idxT;
            } // end yR
        } // end xR
    } // end level
    
    // add offset to motion vectors
    for(size_t i=0; i<matching_table_pyramid[levels].numel(); i++) {
        if(matching_table_pyramid[levels](i)!=MAX_SIZE_T) {
            matching_table_pyramid[levels](i) += offset;
        }
    }
    
    return matching_table_pyramid[levels];
}

/* this function returns the motion field as a volume having size equal to
   the intput video, where each position (yR,xR,tR) in the field contains
   the linearized coordinate of the block most similar to the one having
   coordinate (yR,xR,tR) in frame tR+1 (if dir=FORWARD_TIME) or tR-1 (if 
   dir=INVERSE_TIME); when no motion is detected the value will be MAX_SIZE_T */
volume<size_t> get_motion_field(const volume<real_t> &video, 
        size_t pM, size_t pN, size_t win, real_t threshold, int dir = FORWARD_TIME, 
        int matching_mode = FULL_SEARCH, int estimation_mode = SINGLE_SCALE) {
    // frame will have have no correspondances (MAX_SIZE_T), respectively
    volume<size_t> motion_field = volume<size_t>(video.size(), MAX_SIZE_T);
    
    // ease following processing
    size_t half_win = win/2;
    
    // for each voxel in the video
    for(size_t k=0; k<video.getP(); k++) {
        // get current reference frame (depending on how we traverse the video)
        size_t tR = (dir==FORWARD_TIME)? k : video.getP()-k-1;
        // check if we are not in the boundary frames
        if((tR==0 && dir==INVERSE_TIME) || (tR==video.getP()-1 && dir==FORWARD_TIME)) {
            continue;
        }
        
        // get target frame: tR+1 if FORWARD, and tR-1 if INVERSE
        size_t t = tR + ((dir==FORWARD_TIME)? 1 : -1);
        size_t offset = t*video.getM()*video.getN();
        matrix<size_t> motion_vectors = get_motion_vectors(
                video.get_plane(tR), video.get_plane(t), offset, pM,pN, 
                half_win,threshold, matching_mode, estimation_mode);
        motion_field.set_plane( motion_vectors, tR );
    }
    return motion_field;
}

/* get trajectory of block (yR,xR,tR) concatenating motion vectors contained
   in motion_field; the trajectory has maximum size equal to h; each position
   t in {0,h-1} in the trajectory will be a linearized index in video */
std::vector<size_t> get_trajectory(const volume<size_t> &motion_field,
        size_t idx, size_t h) {
    #if SAFE_CHECK
    if(idx>=motion_field.numel()) {
        throw std::out_of_range("get_trajectory : invalid index");
    }
    #endif
    std::vector<size_t> trajectory;
    trajectory.reserve(h);
    // push coordinate of reference block
    trajectory.push_back( idx );
    for(size_t i=1; i<MIN(h,motion_field.getP()); i++) {
        // get coordinate pointed by current block
        idx = motion_field(idx);
        if(idx>=motion_field.numel()) {
            break;
        }
        // push coordinate if index does not correspond to no motion
        trajectory.push_back( idx );
    }
    return trajectory;
}
std::vector<size_t> get_trajectory(const volume<size_t> &motion_field,
        size_t yR, size_t xR, size_t tR, size_t h) {
    #if SAFE_CHECK
    if(yR>=motion_field.getM() || xR>=motion_field.getN() || 
            tR>=motion_field.getP()) {
        throw std::out_of_range("get_trajectory : invalid indices");
    }
    #endif
    size_t idx = motion_field.get_ind(yR,xR,tR);
    return get_trajectory(motion_field, idx, h);
}


/* interpolate trajectory */
std::vector<std::vector<real_t> > interpolate_trajectory(const volume<size_t> &motion_field, 
        size_t yR, size_t xR, size_t tR, size_t pM, size_t pN, size_t pT, 
        real_t factor, int mode = LINEAR) {
    // extract reference trajectory
    std::vector<size_t> trajectory_ref = get_trajectory(motion_field, yR, xR, tR, pT);
    #if SAFE_CHECK
    if(trajectory_ref.size()!=pT) {
        throw std::logic_error("interpolate_trajectory : invalid trajectory length");
    }
    if(pT==1) {
        throw std::invalid_argument("interpolate_trajectory : invalid trajectory length");
    }
    #endif
    // extract coordinates of reference trajectory
    std::vector<std::vector<real_t> > traj_ind(2);
    std::vector<std::vector<real_t> > traj_ind_intp(5);
    for(size_t n=0; n<pT; n++) {
        traj_ind[0].push_back( motion_field.get_subM(trajectory_ref[n]) );
        traj_ind[1].push_back( motion_field.get_subN(trajectory_ref[n]) );
    }
    // indices of the trajectory in original grid (x)
    traj_ind_intp[4] = range<real_t>(tR,tR+pT-1);
    for(size_t n=std::ceil(traj_ind_intp[4].at(0)*factor); n<=std::floor(traj_ind_intp[4].back()*factor); n++) {
        // indices of trajectory in factor-scaled grid (traj_t_intp)
        traj_ind_intp[2].push_back( n );
        // coordinates of interpolation points in original grid (xi)
        traj_ind_intp[3].push_back( n/factor );
    }
    // interpolate missing positions in reference trajectory
    if(traj_ind_intp[3].size()>0) {
        if(mode==CUBIC && traj_ind[0].size()>2) {
            traj_ind_intp[0] = interp1_pchip(traj_ind_intp[4], traj_ind[0], traj_ind_intp[3]);
            traj_ind_intp[1] = interp1_pchip(traj_ind_intp[4], traj_ind[1], traj_ind_intp[3]);
        } else {
            traj_ind_intp[0] = interp1_linear(traj_ind_intp[4], traj_ind[0], traj_ind_intp[3]);
            traj_ind_intp[1] = interp1_linear(traj_ind_intp[4], traj_ind[1], traj_ind_intp[3]);
        }
        // check boundary conditions and round
        for(size_t n=0; n<traj_ind_intp[0].size(); n++) {
            traj_ind_intp[0][n] = MAX(0,MIN(motion_field.getM()-pM+1,std::round(traj_ind_intp[0][n])));
            traj_ind_intp[1][n] = MAX(0,MIN(motion_field.getN()-pN+1,std::round(traj_ind_intp[1][n])));
        }
    }
    return traj_ind_intp;
    
}

#endif

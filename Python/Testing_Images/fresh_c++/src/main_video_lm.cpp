
#include "mex.h"

#include "master.h"
#include "upsampler_video.hpp"

#if USE_DOUBLE
#define MEX_CLASS mxDOUBLE_CLASS
#else
#define MEX_CLASS mxSINGLE_CLASS
#endif

void mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]) {
    
    // check if input has expected type
    if(USE_DOUBLE && mxIsSingle(prhs[0])) {
        throw std::invalid_argument(" : invalid input class (expected double)");
    }
    if(!USE_DOUBLE && mxIsDouble(prhs[0])) {
        throw std::invalid_argument(" : invalid input class (expected single)");
    }
    if(nrhs!=4) {
        throw std::invalid_argument(" : invalid argument number");
    }
    // check if input has expected number of dimensions
    mwSize ndim = mxGetNumberOfDimensions(prhs[0]);
    if(ndim!=3) {
        throw std::invalid_argument(" : invalid input size (expected 3D grayscale video)");
    }
    // read inputs
    real_t *data_array = (real_t *)mxGetData(prhs[0]);
    std::vector<real_t> data( data_array, data_array+(size_t)mxGetNumberOfElements(prhs[0]) );
    const mwSize *size_data = mxGetDimensions(prhs[0]);
    real_t *data_array_hi = (real_t *)mxGetData(prhs[1]);
    std::vector<real_t> data_hi( data_array_hi, data_array_hi+(size_t)mxGetNumberOfElements(prhs[1]) );
    const mwSize *size_data_hi = mxGetDimensions(prhs[1]);
    
    real_t factorI = (real_t)mxGetScalar(prhs[2]);
    int lm_mode = (int)mxGetScalar(prhs[3]);
    
    // store input data in volume object
    volume<real_t> video = volume<real_t>(size_data[0],size_data[1],size_data[2], data );
    volume<real_t> video_hi = volume<real_t>(size_data_hi[0],size_data_hi[1],size_data_hi[2], data_hi );
    
    // parameters
    size_t pM = 7;
    size_t pN = 7;
    size_t pT = 3;
    
    size_t win_mot = 15;
    real_t threshold_mot = 0.15;
    
    int matching_mode = FULL_SEARCH;
    int estimation_mode = MULTI_SCALE;
    
    size_t K_reg = 8;
    size_t win_reg = 9;
    real_t threshold_reg = 0.3;
    
    size_t pM_lm = 7;
    size_t pN_lm = 7;
    size_t pT_lm = 3;
    size_t K_lm = 4;
    size_t win_lm = 9;
    real_t threshold_lm = 0.3;
    
    // estimate motion fields
    #if VERBOSE
    std::cout << "Motion Estimation" << std::endl;
    #endif
    volume<size_t> MF = get_motion_field(video, pM,pN, win_mot,threshold_mot, 
            FORWARD_TIME, matching_mode, estimation_mode);
    volume<size_t> MF_hi = get_motion_field(video_hi, pM,pN, win_mot,threshold_mot, 
            FORWARD_TIME, matching_mode, estimation_mode);
    
    // create intermediate levels
    #if VERBOSE
    std::cout << "Create Intermediate Levels with Interpolation and Registration" << std::endl;
    #endif
    volume<real_t> video_low, video_inter_intp;
    volume<real_t> video_inter;
    volume<real_t> video_inter_est;
    volume<size_t> MF_inter;
    if(video.getP()==video_hi.getP()) {
        // downsampling to intermediate level
        video_inter = interpolation_along_motion_field(MF, video, pM, pN, pT, factorI/2, LINEAR, 1);
        MF_inter = get_motion_field(video_inter, pM,pN, win_mot,threshold_mot, 
                FORWARD_TIME, matching_mode, estimation_mode);
        video_low = interpolation_along_motion_field(MF_inter, video_inter, pM, pN, pT, 0.5, LINEAR);
        volume<size_t> MF_low = get_motion_field(video_low, pM,pN, win_mot,threshold_mot, 
                FORWARD_TIME, matching_mode, estimation_mode);
        video_inter_intp = interpolation_along_motion_field(MF_low, video_low, pM, pN, pT, 2, LINEAR);
        video_inter_intp.set_plane(video_inter.get_plane(video_inter.getP()-1),video_inter_intp.getP()-1);
        // upsampling by 2 with nl registration
        video_inter_est = upsample_by_two_using_nl_registration(MF_low, video_inter_intp, video_low, 
                K_reg, pM, pN, pT, win_reg, threshold_reg);
    } else {
        // upsampling to intermediate level by interpolation along motion field
        video_inter = interpolation_along_motion_field(MF, video, pM, pN, pT, factorI, LINEAR, 1);
        MF_inter = get_motion_field(video_inter, pM,pN, win_mot,threshold_mot, 
                FORWARD_TIME, matching_mode, estimation_mode);
        video_low = interpolation_along_motion_field(MF_inter, video_inter, pM, pN, pT, 0.5, LINEAR);
        volume<size_t> MF_low = get_motion_field(video_low, pM,pN, win_mot,threshold_mot, 
                FORWARD_TIME, matching_mode, estimation_mode);
        video_inter_intp = interpolation_along_motion_field(MF_low, video_low, pM, pN, pT, 2, LINEAR);
        video_inter_intp.set_plane(video_inter.get_plane(video_inter.getP()-1),video_inter_intp.getP()-1);
        video_inter_est = upsample_by_two_using_nl_registration(MF_low, video_inter_intp, video_low, 
                K_reg, pM, pN, pT, win_reg, threshold_reg);
    }
    
    // lm
    volume<real_t> video_lm;
    if(lm_mode) {
        #if VERBOSE
        std::cout << "Linear Mapping with Motion Estimation" << std::endl;
        #endif
        video_lm = correction_by_linear_mapping_with_motion_field(
                MF_hi, MF_inter, video_hi, video_inter, video_inter_est,
                K_lm, pM, pN, pT, win_lm, threshold_lm);
    } else {
        #if VERBOSE
        std::cout << "Linear Mapping without Motion Estimation" << std::endl;
        #endif
        video_lm = correction_by_linear_mapping(
                video_hi, video_inter, video_inter_est, 
                K_lm, pM_lm, pN_lm, pT_lm, win_lm, threshold_lm);
    }
    
    // write output
    mwSize *nsize = new mwSize[ndim];
    nsize[0] = video_lm.getM();
    nsize[1] = video_lm.getN();
    nsize[2] = video_lm.getP();
    plhs[0] = mxCreateNumericArray(ndim, nsize, MEX_CLASS, mxREAL);
    real_t *data_out_lm_res = (real_t *)mxGetData(plhs[0]);
    for(size_t i=0; i<video_lm.numel(); i++) {
        data_out_lm_res[i] = video_lm(i);
    }
    
    nsize[0] = video_inter.getM();
    nsize[1] = video_inter.getN();
    nsize[2] = video_inter.getP();
    plhs[1] = mxCreateNumericArray(ndim, nsize, MEX_CLASS, mxREAL);
    real_t *data_out_inter = (real_t *)mxGetData(plhs[1]);
    for(size_t i=0; i<video_inter.numel(); i++) {
        data_out_inter[i] = video_inter(i);
    }
    
    nsize[0] = video_inter_est.getM();
    nsize[1] = video_inter_est.getN();
    nsize[2] = video_inter_est.getP();
    plhs[2] = mxCreateNumericArray(ndim, nsize, MEX_CLASS, mxREAL);
    real_t *data_out_inter_est = (real_t *)mxGetData(plhs[2]);
    for(size_t i=0; i<video_inter_est.numel(); i++) {
        data_out_inter_est[i] = video_inter_est(i);
    }
    
    nsize[0] = video_inter_intp.getM();
    nsize[1] = video_inter_intp.getN();
    nsize[2] = video_inter_intp.getP();
    plhs[3] = mxCreateNumericArray(ndim, nsize, MEX_CLASS, mxREAL);
    real_t *data_out_low = (real_t *)mxGetData(plhs[3]);
    for(size_t i=0; i<video_inter_intp.numel(); i++) {
        data_out_low[i] = video_inter_intp(i);
    }
    
    free(nsize);
    
}
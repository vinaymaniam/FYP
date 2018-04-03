
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
    if(nrhs!=2) {
        throw std::invalid_argument(" : invalid argument number");
    }
    // check if input has expected number of dimensions
    mwSize ndim = mxGetNumberOfDimensions(prhs[0]);
    const mwSize *size_data = mxGetDimensions(prhs[0]);
    if(ndim!=3) {
        throw std::invalid_argument(" : invalid input size (expected grayscale video)");
    }
    // read inputs
    real_t *data_array = (real_t *)mxGetData(prhs[0]);
    std::vector<real_t> data( data_array, data_array+(size_t)mxGetNumberOfElements(prhs[0]) );
    real_t factorT = (real_t)mxGetScalar(prhs[1]);
    
    // store input data in volume object
    size_t M = size_data[0];
    size_t N = size_data[1];
    size_t T = size_data[2];
    volume<real_t> video = volume<real_t>(M,N,T, data );
    
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
    
    // motion field
    #if VERBOSE
    std::cout << "Motion Estimation" << std::endl;
    #endif
    // interpolation along motion field
    volume<size_t> MF = get_motion_field(video, pM,pN, win_mot,threshold_mot, 
            FORWARD_TIME, matching_mode, estimation_mode);
    #if VERBOSE
    std::cout << "Interpolation along Motion Field" << std::endl;
    #endif
    volume<real_t> video_intp = interpolation_along_motion_field(MF, video, pM, pN, pT, factorT, LINEAR);
    
    
    // upsampled video via nonlocal patch registration
    #if VERBOSE
    std::cout << "Patch Registration with Motion Estimation" << std::endl;
    #endif
    volume<real_t> video_nl = upsample_by_two_using_nl_registration(MF, video_intp, video,
            K_reg, pM, pN, pT, win_reg, threshold_reg);
    
    // output 0 -> motion field
    mwSize *nsize = new mwSize[ndim];
    nsize[0] = MF.getM();
    nsize[1] = MF.getN();
    nsize[2] = MF.getP();
    plhs[0] = mxCreateNumericArray(ndim, nsize, mxUINT64_CLASS, mxREAL);
    uint64_T *data_out_mf = (uint64_T *)mxGetData(plhs[0]);
    for(size_t i=0; i<MF.numel(); i++) {
        data_out_mf[i] = (uint64_T)MF(i);
    }
    
    // output 1 -> registered input
    //mwSize *nsize = new mwSize[ndim];
    nsize[0] = video_nl.getM();
    nsize[1] = video_nl.getN();
    nsize[2] = video_nl.getP();
    plhs[1] = mxCreateNumericArray(ndim, nsize, MEX_CLASS, mxREAL);
    real_t *data_out_nl_res = (real_t *)mxGetData(plhs[1]);
    for(size_t i=0; i<video_nl.numel(); i++) {
        data_out_nl_res[i] = video_nl(i);
    }
    
    // output 2 -> interpolated input
    //mwSize *nsize = new mwSize[ndim];
    nsize[0] = video_intp.getM();
    nsize[1] = video_intp.getN();
    nsize[2] = video_intp.getP();
    plhs[2] = mxCreateNumericArray(ndim, nsize, MEX_CLASS, mxREAL);
    real_t *data_out_intp = (real_t *)mxGetData(plhs[2]);
    for(size_t i=0; i<video_intp.numel(); i++) {
        data_out_intp[i] = video_intp(i);
    }
    free(nsize);
    
}

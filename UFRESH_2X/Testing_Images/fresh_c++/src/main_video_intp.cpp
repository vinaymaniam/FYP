
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
    if(!mxIsClass(prhs[0],"uint64")) {
        throw std::invalid_argument(" : invalid input class (expected uint64)");
    }
    if(USE_DOUBLE && mxIsSingle(prhs[1])) {
        throw std::invalid_argument(" : invalid input class (expected double)");
    }
    if(!USE_DOUBLE && mxIsDouble(prhs[1])) {
        throw std::invalid_argument(" : invalid input class (expected single)");
    }
    if(nrhs!=3) {
        throw std::invalid_argument(" : invalid argument number");
    }
    // check if input has expected number of dimensions
    mwSize ndim = mxGetNumberOfDimensions(prhs[0]);
    const mwSize *size_data = mxGetDimensions(prhs[0]);
    if(ndim!=3) {
        throw std::invalid_argument(" : invalid input size (expected grayscale video)");
    }
    // read inputs
    uint64_T *data_array_mf = (uint64_T *)mxGetData(prhs[0]);
    std::vector<uint64_T> data_mf( data_array_mf, data_array_mf+(size_t)mxGetNumberOfElements(prhs[0]) );
    real_t *data_array = (real_t *)mxGetData(prhs[1]);
    std::vector<real_t> data( data_array, data_array+(size_t)mxGetNumberOfElements(prhs[1]) );
    real_t factorT = (real_t)mxGetScalar(prhs[2]);
    
    #if SAFE_CHECK
    if(data_mf.size()!=data.size()) {
        throw std::invalid_argument(" : size of video and motion field do not agree");
    }
    #endif
    
    // store input data in volume object
    size_t M = size_data[0];
    size_t N = size_data[1];
    size_t T = size_data[2];
    volume<real_t> video = volume<real_t>(M,N,T, data );
    volume<size_t> MF = volume<size_t>(M,N,T);
    for(size_t i=0; i<M*N*T; i++) {
        MF(i) = (uint64_T)data_mf.at(i);
    }
    
    // parameters
    size_t pM = 7;
    size_t pN = 7;
    size_t pT = 3;
    
    #if VERBOSE
    std::cout << "Interpolation along Motion Field" << std::endl;
    #endif
    volume<real_t> video_intp = interpolation_along_motion_field(MF, video, pM, pN, pT, factorT, LINEAR);
    
    // output 0 -> interpolated input
    mwSize *nsize = new mwSize[ndim];
    nsize[0] = video_intp.getM();
    nsize[1] = video_intp.getN();
    nsize[2] = video_intp.getP();
    plhs[0] = mxCreateNumericArray(ndim, nsize, MEX_CLASS, mxREAL);
    real_t *data_out_intp = (real_t *)mxGetData(plhs[0]);
    for(size_t i=0; i<video_intp.numel(); i++) {
        data_out_intp[i] = video_intp(i);
    }
    free(nsize);
    
}

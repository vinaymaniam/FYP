
// compile as stand-alone (reading file from disk) or as Matlab mex interface
#ifndef MEX_COMPILE_FLAG
#define MEX_COMPILE_FLAG 0
#endif

#include <stdexcept>
#include <vector>

#include "master.h"
#include "wavelet.hpp"
#include "upsampler.hpp"
#include "io.hpp"

#if MEX_COMPILE_FLAG

#include "mex.h"
#if USE_DOUBLE
#define MEX_CLASS mxDOUBLE_CLASS
#else
#define MEX_CLASS mxSINGLE_CLASS
#endif

void mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]) {
#else

int main(int argc, char *argv[]) {
#endif
    
    #if MEX_COMPILE_FLAG
    // read input variables from Matlab call
    if(USE_DOUBLE && mxIsSingle(prhs[0])) {
        mexErrMsgTxt("Invalid input class (expected double)");
    }
    if(!USE_DOUBLE && mxIsDouble(prhs[0])) {
        mexErrMsgTxt("Invalid input class (expected single)");
    }
    if(nrhs!=6) {
        mexErrMsgTxt("Invalid argument number");
    }
    // read from matlab input
    real_t *data_array = (real_t *)mxGetData(prhs[0]);
    std::vector<real_t> data( data_array, data_array+(size_t)mxGetNumberOfElements(prhs[0]) );
    mwSize ndim = mxGetNumberOfDimensions(prhs[0]);
    const mwSize *size_data = mxGetDimensions(prhs[0]);
    size_t M = size_data[0];
    size_t N = size_data[1];
    size_t channels = (ndim==3)? 3 : 1;
    std::vector<matrix<real_t> > I_rgb;
    I_rgb.reserve(channels);
    for(size_t i=0; i<channels; i++) {
        std::vector<real_t> data_c(data.begin()+M*N*i,data.begin()+M*N*(i+1));
        I_rgb.push_back( matrix<real_t>(M,N,data_c) );
    }
    size_t level        = (size_t)mxGetScalar(prhs[1]);
    int upsampling_mode = (int)mxGetScalar(prhs[2]);
    int correction_mode = (int)mxGetScalar(prhs[3]);
    int algorithm_mode  = (int)mxGetScalar(prhs[4]);
    int profile         = (int)mxGetScalar(prhs[5]);
    #else
    if(argc!=8) {
        std::cout << 
                "FRESH usage:" << std::endl <<
                "./fresh input output level diag_reg lin_map test_mode fast" << std::endl <<
                "       input      :  path to the input image (only BMP is supported) - string" << std::endl <<
                "       output     :  path of the output image (only BMP is supported) - string" << std::endl <<
                "       level      :  upsampling factor is 2^level - integer {1,2,...}" << std::endl <<
                "       diag_reg   :  enable diagonal regularization - boolean {0,1}" << std::endl <<
                "       lin_map    :  enable linear mapping - boolean {0,1}" << std::endl <<
                "       test_mode  :  downsample input before upsampling - boolean {0,1}" << std::endl <<
                "       fast       :  enable fast profile - boolean {0,1}" << std::endl <<
                "The upsampled image will be written to the filepath specified by \"output\"." << std::endl;
        return 1;
    }
    // read parameters
    std::string in_filename = argv[1];
    std::string out_filename = argv[2];
    size_t level        = (size_t)atoi(argv[3]);
    int upsampling_mode = (int)atoi(argv[4]);
    int correction_mode = (int)atoi(argv[5]);
    int algorithm_mode  = (int)atoi(argv[6]);
    int profile         = (int)atoi(argv[7]);
    
    // read image from file
    image_t<real_t> I_rgb;
    try {
        I_rgb = imread(in_filename);
    } catch(std::exception &e) {
        std::cout << std::string("An error occurred while reading input image. \n\nFRESH::") + 
                e.what() << std::endl;
        return 1;
    }
    #endif
    
    // FRI filter parameters
    int filter = RBIO2_8;
    int filter_fri = SPLN4;
    
    // upsampling scaling factor: M-by-N -> (factor*M)-by-(factor*N)
    // being factor = 2^level, where level is an input argument
    real_t factor = (real_t)(1<<level);

    // output variable for upsampled data
    image_t<real_t> I_rgb_up;
    try {
        // high-resolution input, thus we first artificially downsample it
        // and then we upsample its low-pass subband (in this way we can
        // compare the original hi-res versions with the upsampled one)
        if(algorithm_mode==HI_RES_IN) {
            filter = BIOR4_4;
            filter_fri = BIOR4_4;
            wavelet<real_t> wlt(filter);
            for(size_t i=0; i<I_rgb.size(); i++) {
                std::vector<matrix<real_t> > c = wlt.wavedec2(I_rgb[i],level);
                I_rgb[i] = wlt.appcoef2(c,level) / factor;
            }
        }
        // convert to YUV and extract luminance
        image_t<real_t> I_yuv = (I_rgb.size()==3)? rgb2ycbcr(I_rgb, ROUND) : I_rgb;
        matrix<real_t> I_luma = I_yuv[0];
        // init singleton upsampler object
        upsampler UPSMPLR(filter, filter_fri, 
                upsampling_mode, correction_mode, algorithm_mode, profile);
        // call upsampler method
        #if PARALLEL
        matrix<real_t> I_up_n = UPSMPLR.upsample_parallel(I_luma, level);
        #else
        matrix<real_t> I_up_n = UPSMPLR.upsample(I_luma, level);
        #endif
        I_rgb_up.reserve(I_rgb.size());
        I_up_n.round();
        I_rgb_up.push_back(I_up_n);
        for(size_t i=1; i<I_rgb.size(); i++) {
            // upsampling chrominance channels (if any) by cubic interpolation
            matrix<real_t> I_upres_crcb = UPSMPLR.shifted_resize(I_yuv[i], factor);
            I_upres_crcb.round();
            I_rgb_up.push_back(I_upres_crcb);
        }
        // convert back to RGB
        if(I_rgb.size()==3) {
            I_rgb_up = ycbcr2rgb(I_rgb_up, ROUND);
        }
        // clip to uint8 range
        for(size_t i=0; i<I_rgb_up.size(); i++) {
            I_rgb_up[i].clip(0, 255);
        }
    } catch(std::exception &e) {
        const std::string err_msg = 
                std::string("An error occurred during upsampling. \n\nFRESH::") + e.what();
        #if MEX_COMPILE_FLAG==1
        mexErrMsgTxt(err_msg.c_str());
        #else
        std::cout << err_msg << std::endl;
        return 1;
        #endif
    }
    
    #if MEX_COMPILE_FLAG
    // link high-res image to Matlab first output variable
    mwSize *nsize = new mwSize[ndim];
    nsize[0] = (mwSize)I_rgb_up[0].getM();
    nsize[1] = (mwSize)I_rgb_up[0].getN();
    if(ndim==3) {
        nsize[2] = 3;
    }
    plhs[0] = mxCreateNumericArray(ndim, nsize, MEX_CLASS, mxREAL);
    real_t *data_out_hi_res = (real_t *)mxGetData(plhs[0]);
    size_t nelem = prod1(I_rgb_up[0].size());
    for(size_t i=0; i<channels; i++) {
        for(size_t j=0; j<nelem; j++) {
            data_out_hi_res[j+i*nelem] = I_rgb_up[i](j);
        }
    }
    // link bicubic-resized image to Matlab second output variable
    image_t<real_t> I_bic;
    I_bic.reserve(I_rgb.size());
    for(size_t i=0; i<I_rgb.size(); i++) {
        I_bic.push_back( resize(I_rgb[i], factor, BICUBIC) );
    }
    nsize[0] = (mwSize)I_bic[0].getM();
    nsize[1] = (mwSize)I_bic[0].getN();
    if(ndim==3) {
        nsize[2] = 3;
    }
    plhs[1] = mxCreateNumericArray(ndim, nsize, MEX_CLASS, mxREAL);
    real_t *data_out_bic_res = (real_t *)mxGetData(plhs[1]);
    nelem = prod1(I_bic[0].size());
    for(size_t i=0; i<channels; i++) {
        for(size_t j=0; j<nelem; j++) {
            data_out_bic_res[j+i*nelem] = I_bic[i](j);
        }
    }
    // link low-res image to Matlab third output variable
    nsize[0] = (mwSize)I_rgb[0].getM();
    nsize[1] = (mwSize)I_rgb[0].getN();
    if(ndim==3) {
        nsize[2] = 3;
    }
    plhs[2] = mxCreateNumericArray(ndim, nsize, MEX_CLASS, mxREAL);
    real_t *data_out_lo_res = (real_t *)mxGetData(plhs[2]);
    nelem = prod1(I_rgb[0].size());
    for(size_t i=0; i<channels; i++) {
        for(size_t j=0; j<nelem; j++) {
            data_out_lo_res[j+i*nelem] = I_rgb[i](j);
        }
    }
    delete [] nsize;
    #else
    try {
        // saving upsampled image to BMP file
        imwrite(I_rgb_up, out_filename);
    } catch(std::exception &e) {
        std::cout << std::string("An error occurred while writing output image. \n\nFRESH::") + 
                e.what() << std::endl;
        return 1;
    }
	return 0;
    #endif
    
}


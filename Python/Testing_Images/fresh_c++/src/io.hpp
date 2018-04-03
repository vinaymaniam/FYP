
#ifndef IO_H
#define IO_H

#include <iostream>
#include <stdexcept>
#include <vector>

#include "master.h"
#include "matrix.hpp"
#include "image.hpp"

#ifndef WITH_JPEG
#define WITH_JPEG 0
#endif

#if WITH_JPEG
extern "C"{
#include <jpeglib.h>
}
#endif

typedef int LONG;
typedef unsigned short WORD;
typedef unsigned int DWORD;
typedef unsigned char BYTE;

#pragma pack(2)
typedef struct {
    WORD  bfType;          // file type (must be BM)
    DWORD bfSize;          // file size (in bytes)
    WORD  bfReserved1;     // reserved (must be 0)
    WORD  bfReserved2;     // reserved (must be 0)
    DWORD bfOffBits;       // offset (in bytes) from header to the bitmap
} BITMAPFILEHEADER;
#pragma pack()

typedef struct {
    DWORD biSize;          // size (in bytes) required by the structure
    LONG  biWidth;         // width of bitmap (in pixels)
    LONG  biHeight;        // height of bitmap (in pixels)
    WORD  biPlanes;        // number of planes for the target device (must be 1)
    WORD  biBitCount;      // bit-per-pixels (2^biBitCount)
    DWORD biCompression;   // compression of the bitmap
    DWORD biSizeImage;     // size of image (in bytes)
    LONG  biXPelsPerMeter; // x-resolution of target device
    LONG  biYPelsPerMeter; // y-resolution of target device
    DWORD biClrUsed;       // color actually used in the bitmap
    DWORD biClrImportant;  // number of colors required to display the image
} BITMAPINFOHEADER;

typedef struct {
    BYTE rgbBlue;
    BYTE rgbGreen;
    BYTE rgbRed;
    BYTE rgbReserved;
} RGBQUAD;

/* read BMP file and return as output arguments a vector containing three 
   matrix object (the channels of the image if the file is RGB), or a
   vector containing a single matrix object (if the image is grayscale) */
std::vector<matrix<real_t> > read_bmp(const std::string &filename) {
    // read file and check if the file exist
    FILE* file = fopen(filename.c_str(), "rb");
    if(!file) {
        throw std::invalid_argument("read_bmp : error reading input file " + filename);
    }
    
    // read headers
    BITMAPFILEHEADER file_header;
    BITMAPINFOHEADER info_header;
    fread(&file_header, sizeof(BITMAPFILEHEADER), 1, file);
    fread(&info_header, sizeof(BITMAPINFOHEADER), 1, file);
    
    // check if the file is an actual BMP file, i.e. bfType should correspond 
    // to the letters 'B' (0×42) and 'M' (0x4D), and no compression is used
    // !!! be aware if the system is big or little endian !!!
	if(file_header.bfType!=0x4D42 || info_header.biCompression!=0) {
		throw std::invalid_argument("read_bmp : invalid BMP file");
	}
    // check if the expected number of bit-per-pixel is met
    if(info_header.biBitCount!=8 && info_header.biBitCount!=24) {
        throw std::invalid_argument("read_bmp : invalid BMP bit-per-pixel value " + 
                std::to_string(info_header.biBitCount));
    }
    size_t N = info_header.biWidth;
    size_t M = info_header.biHeight;
    size_t C = (info_header.biBitCount)>>3; // color channels (either 1 or 3)
    // width of bitmap N must be divisible by 4, if not a padding is applied
    size_t P = 4 - ((N * C) % 4);
    P = (P==4)? 0 : P;
    
    // read bitmap and store to buffer
    unsigned char *buffer = new unsigned char[(N*C+P)*M];
    fseek(file, file_header.bfOffBits, SEEK_SET);
	fread(buffer, sizeof(unsigned char), (N*C+P)*M, file);
    fclose(file);
    
    // bitmap is row-major bottom-up from bottom-left pixel, thus width of 
    // bitmap is height in matrix, and viceveresa
    // create output matrix
    std::vector<matrix<real_t> > img;
    img.reserve(C);
    for(size_t k=0; k<C; k++) {
        img.push_back( matrix<real_t>(M,N) );
    }
    for(size_t i=0; i<M; i++) {
        for(size_t j=0; j<N; j++) {
            for(size_t k=0; k<C; k++) {
                // each pixel is a single byte (for grayscale bitmaps), or
                // a sequence of three bytes (for RGB images) which
                // correspond to blue green and red channel
                int val = (int)buffer[C*(i*N+j)+i*P+k];
                img[C-k-1](M-i-1,j) = (real_t)val;
            }
        }
    }
    delete [] buffer;
    return img;
}

/* write image to BMP file, the input image is a vector containing three 
   matrix object (the channels of the image) */
void write_bmp(const std::vector<matrix<real_t> > &img, const std::string &filename) {
    // getting image size 
    size_t M = img[0].getM();
    size_t N = img[0].getN();
    size_t C = img.size();
    // width of bitmap M must be divisible by 4, if not a padding is applied
    size_t P = 4 - ((N * C) % 4);
    P = (P==4)? 0 : P;
    
    // define size of headers and bitmap
    size_t headers_size = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER) + 
            (C==1)*sizeof(RGBQUAD)*256; // add colormap if grayscale
    size_t image_size = (N*C+P)*M;
    
    // define headers
    BITMAPFILEHEADER file_header = {0};
    file_header.bfType = (WORD)0x4D42; //'B' + ('M'<<8);
    file_header.bfOffBits = (DWORD)headers_size;
    file_header.bfSize = (DWORD)(headers_size + image_size);
    
    BITMAPINFOHEADER info_header = {0};
    info_header.biSize = (DWORD)sizeof(BITMAPINFOHEADER);
    info_header.biBitCount = (WORD)(C*8);
    info_header.biClrImportant = 0;
    info_header.biClrUsed = (DWORD)(C==3? 0 : 256);
    info_header.biCompression = 0;
    info_header.biHeight = (LONG)M;
    info_header.biWidth = (LONG)N;
    info_header.biPlanes = (WORD)1;
    info_header.biSizeImage = (DWORD)image_size;
    
    // open and link output file
    FILE * outfile;
    if((outfile = fopen(filename.c_str(), "wb"))==NULL) {
        throw std::invalid_argument("write_bmp : error writing output file " + filename);
    }
    
    // write headers
    fwrite(&file_header, sizeof(BITMAPFILEHEADER), 1, outfile);
    fwrite(&info_header, sizeof(BITMAPINFOHEADER), 1, outfile);
    
    // if grayscale we should add the colormap before the data
    if(C==1) {
        RGBQUAD color_map[256];
        for(size_t i=0; i<256; i++) {
            color_map[i].rgbBlue = (BYTE)i;
            color_map[i].rgbGreen = (BYTE)i;
            color_map[i].rgbRed = (BYTE)i;
            color_map[i].rgbReserved = (BYTE)0;
        }
        fwrite(&color_map, sizeof(RGBQUAD)*256, 1, outfile);
    }
    
    // write data to file (row-major bottom-up) in BGR form
    for(size_t i=0; i<M; i++) {
        for(size_t j=0; j<N; j++) {
            for(size_t k=0; k<C; k++) {
                unsigned char val = (unsigned char)((int)img[C-k-1](M-i-1,j));
                fwrite(&val, sizeof(unsigned char), 1, outfile);
            }
        }
        if(P>0) {
            unsigned char zero = (unsigned char)0;
            for(size_t n=0; n<P; n++) {
                fwrite(&zero, sizeof(unsigned char), 1, outfile);
            }
        }
    }
    fclose(outfile);
}

#if WITH_JPEG
/* custom bind to error function called when jpeglib returns error */
void jpeg_error_exit(j_common_ptr cinfo) {
    char jpegLastErrorMsg[JMSG_LENGTH_MAX];
    // create the message 
    ( *( cinfo->err->format_message ) ) ( cinfo, jpegLastErrorMsg );
    // jump to the setjmp point 
    throw std::runtime_error( jpegLastErrorMsg );
}

/* read JPEG file and return as output arguments a vector containing three 
   matrix object (the channels of the image); this function uses libjpeg 
   NOTE: in order to replicate Matlab imread, libjpeg v6b is required */
std::vector<matrix<real_t> > read_jpeg(const std::string &filename) {
    struct jpeg_decompress_struct info;
    struct jpeg_error_mgr err;
    
    FILE* file = fopen(filename.c_str(), "rb");
    
    if(!file) {
        throw std::invalid_argument("read_jpeg : error reading input file " + filename);
    }
    
    info.err = jpeg_std_error(&err);
    err.error_exit = jpeg_error_exit;
    jpeg_create_decompress(&info);
    jpeg_stdio_src(&info, file);
    jpeg_read_header(&info, TRUE);
    jpeg_start_decompress(&info);
    
    //set width and height, channels, and data_size
    size_t M = info.output_height;
    size_t N = info.output_width;
    size_t channels = (size_t)info.num_components;
    
    // output arguments
    std::vector<matrix<real_t> > img;
    img.reserve(channels);
    for(size_t c=0; c<channels; c++) {
        img.push_back( matrix<real_t>(M, N) );
    }
    
    // row is an array containing the current scanned rows in each channel
    // of the image; the values of the channels are interleaved as
    // {R1, G1, B1, R2, G2, B2, ...., RN GN BN}, N being the width of the image
    unsigned char *row = new unsigned char[N*channels];
    matrix<real_t> scan = matrix<real_t>(N, channels);
    while(info.output_scanline < M) {
        jpeg_read_scanlines(&info, &row, 1);
        // read values in scanline for all channels
        for(size_t c=0; c<channels; c++) {
            for(size_t k=0; k<N; k++) {
                // convert unsigend char to int and then to real_t
                scan(k,c) = (real_t)((int)row[3*k+c]);
            }
            // save current row to proper channel in output arg
            img[c].set_row(scan.get_col(c), info.output_scanline-1);
        }
    }
    
    // release resources
    fclose(file);
    jpeg_finish_decompress(&info);
    delete [] row;
    
    return img;
}

/* write image to JPEG file, the input image is a vector containing three 
   matrix object (the channels of the image) */
void write_jpeg(const std::vector<matrix<real_t> > &img, 
        const std::string &filename, int quality = 100) {
    // libjpeg structures
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;
    
    // set up error handler and JPEG compression object
    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);
    
    // open and link output file
    FILE * outfile;
    if((outfile = fopen(filename.c_str(), "wb"))==NULL) {
        throw std::invalid_argument("write_jpeg : error writing output file " + filename);
    }
    jpeg_stdio_dest(&cinfo, outfile);
    
    // compression parameters
    cinfo.image_width = (unsigned int)img[0].getN();
    cinfo.image_height = (unsigned int)img[0].getM();
    cinfo.input_components = (int)img.size(); // number of channels
    cinfo.in_color_space = JCS_RGB; // colorspace
    // other parameters set as default
    jpeg_set_defaults(&cinfo);
    // set desired JPG compression quality (TRUE limits to baseline-JPEG,
    // i.e. the computed quantization table entries are clipped in [1,255])
    // quality can be in range [0,100]
    jpeg_set_quality(&cinfo, quality, TRUE);
    
    // start compression (TRUE ensures that a complete interchange-JPEG
    // datastream will be written and all Huffman tables shall be emited)
    jpeg_start_compress(&cinfo, TRUE);
    
    // row is an array containing a row in each channel of the image; the 
    // values are interleaved as {R1, G1, B1, R2, G2, B2, ...., RN GN BN},
    // N being the width (cinfo.image_width) of the image
    unsigned char *row = new unsigned char[cinfo.image_width*cinfo.input_components];
    for(size_t i=0; i<cinfo.image_height; i++) {
        for(size_t j=0; j<cinfo.image_width; j++) {
            for(size_t c=0; c<cinfo.input_components; c++) {
                // convert to int and then to unsigend char
                row[3*j+c] = (unsigned char)((int) img[c](i,j));
            }
        }
        // write current line to output
        jpeg_write_scanlines(&cinfo, &row, 1);
    }
    
    // release resources
    jpeg_finish_compress(&cinfo);
    fclose(outfile);
    jpeg_destroy_compress(&cinfo);
}
#endif

/* read any supported image file (currently either JPEG or BMP) to a vector 
   containing three matrix objects (if the image is color), or a vector 
   containing a single matrix object (if the image is grayscale) */
std::vector<matrix<real_t> > imread(const std::string &filename) {
    std::string ext = filename.substr(filename.find_last_of(".")+1);
    if(ext=="bmp") {
        return read_bmp(filename);
    #if WITH_JPEG
    } else if(ext=="jpeg" || ext=="jpg") {
        return read_jpeg(filename);
    #endif
    } else {
        throw std::invalid_argument("imread : unsupported file type " + ext);
    }
}

/* write vector of matrices to file; the element in the vectors are the 
   channels of the image and thus can be either 1 (grayscale data) or 3 
   (color data); output file is determined by the extension of filename
   and can be any supported image format (currently either JPEG or BMP) */
void imwrite(const std::vector<matrix<real_t> > &img, const std::string &filename) {
    #if SAFE_CHECK
    if(img.size()!=1 && img.size()!=3) {
        throw std::invalid_argument("imwrite : invalid input dimensions");
    }
    #endif
    std::string ext = filename.substr(filename.find_last_of(".")+1);
    if(ext=="bmp") {
        return write_bmp(img, filename);
    #if WITH_JPEG
    } else if(ext=="jpeg" || ext=="jpg") {
        return write_jpeg(img, filename);
    #endif
    } else {
        throw std::invalid_argument("imwrite : unsupported file type " + ext);
    }
}

#endif


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                COMPILING                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
with_jpeg = 0;
parallel = 0;
libjpeg_dir = '/usr/local/Cellar/jpeg/9a/'; % path to libjpeg directory
if with_jpeg && exist(libjpeg_dir,'dir')
    libjpeg = ['-I',libjpeg_dir,'include/ -L',libjpeg_dir,'lib/ -ljpeg -DWITH_JPEG=1 '];
else
    libjpeg = '-DWITH_JPEG=0 ';
end
if isunix
    coptimflags = '-O3 -ffast-math';
else
    coptimflags = '/O2gtx /fp:fast';
end
eval(['mex COPTIMFLAGS="$COPTIMFLAGS ',coptimflags,'" ',...
    ' -DMEX_COMPILE_FLAG=1 -DPARALLEL=',int2str(parallel),' ',...
    libjpeg,' -output fresh src/main.cpp -lmwlapack -lmwblas' ]);
pcode('fresh_gui.m');
clear with_jpeg libjpeg_dir libjpeg coptimflags

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  NOTICE                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Demo software for FRESH Single-Image Super-Resolution
% 
% The software implements the super-resolution algorithm described in:
% [1] X. Wei and P. L. Dragotti, "Sampling piecewise smooth signals and its
%     application to image up-sampling", in ICIP, Sep. 2015, pp. 4293?4297. 
% [2] X. Wei and P. L. Dragotti, "FRESH - FRI-based single-image 
%     super-resolution algorithm", in IEEE Transactions on Image Processing, 
%     vol. 25, no. 8, pp. 3723-3735, Aug. 2016.
% 
% Authors:      Xiaoyao Wei           ivy.wei@imperial.ac.uk
%               Matteo Maggioni       m.maggioni@imperial.ac.uk
%               Pier Luigi Dragotti   p.dragotti@imperial.ac.uk
% 
% Copyright (c) 2015-2016 Imperial College London.
% All rights reserved.
% 
% Disclaimer:
% This work should be used for nonprofit purposes only. Any unauthorized 
% use of these routines for industrial or profit-oriented activities is 
% expressively prohibited. By downloading and/or using any of these files, 
% you implicitly agree to all the terms of the limited license included in 
% this package (Legal_Notice.txt).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             DEMO PARAMETERS                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

directory = 'data'; 
pattern = '*.bmp';

XpathCell = glob(directory, pattern);
Xcell = load_images( XpathCell );

%filename     = 'training/0170x4.png';
test_mode    = 1; % enable or disable test mode:
                  % if test_mode is 1, then we fist downsample the input
                  %   image referenced by `filename', and then we apply the
                  %   desired upsampling `level' to eventually produce an 
                  %   output image having the same size as the input, this
                  %   allows to objectively asses the performance of the method;       
                  % if test_mode is 0, then the upsampled result will have 
                  %   dimension equal to size(input)*2^level, and thus in   
                  %   this case no objective comparison is possible
save_to_file = 1; % save upsampled result to PNG file

level        = 2; % upsample input by scaling factor 2^level (e.g., set 
                  % level=1 to double the size of the image)
diag_reg     = 0; % enable diagonal regularization
lin_map      = 0; % enable error correction by linear mapping
profile      = 0; % computational profile {0->standard, 1->fast}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     DO NOT MODIFY BELOW THIS POINT                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n=1:length(Xcell)
img=Xcell{n};

% info = imfinfo(filename);
% bitdepth = info.BitDepth;
% if strcmpi(info.ColorType,'truecolor')
%     bitdepth = bitdepth/3;
% end
%img = double(imread(filename))/(2^bitdepth-1)*255;
%img=imresize(I,1/2^level);
%fprintf('FRI %dX UPSAMPLING - ', 2^level)
% if test_mode
%     fprintf('%dx%d -> %dx%d \n',...
%         ceil(size(img,1)/2^level),ceil(size(img,2)/2^level), size(img,1),size(img,2))
% else
%     fprintf('%dx%d -> %dx%d \n', ...
%         size(img,1),size(img,2), size(img,1)*2^level,size(img,2)*2^level)
% end

% output parameters
%  img_fresh : the upsampling result obtained via the proposed method
%  img_bic   : the upsampling result using bicubic interpolation (for comparison)
%  img_low   : the low-resolution image (which is equal to the input image 
%              if test_mode is 0)
tic;
[img_fresh, img_bic, img_low] = fresh( img, level, diag_reg, lin_map, test_mode, profile );
toc

%%
% show results
figure,
subplot(2,2,1),imshow(uint8(img)),title('ORIGINAL')
subplot(2,2,2),imshow(uint8(img_low)),title('LOW-RES')
subplot(2,2,3),imshow(uint8(img_bic))
if test_mode
    psnr_bic = 10*log10(255^2/mean((img_bic(:)-img(:)).^2));
    ssim_bic = ssim(img_bic, img);
    title(sprintf('BICUBIC %dX \n PSNR: %.2fdB, SSIM: %.4f',...
        2^level, psnr_bic, ssim_bic))
else
    title('BICUBIC')
end
subplot(2,2,4),imshow(uint8(img_fresh))
if test_mode
    psnr_fresh = 10*log10(255^2/mean((img_fresh(:)-img(:)).^2));
    ssim_fresh = ssim(img_fresh, img);
    title(sprintf('FRESH %dX \n PSNR: %.2fdB, SSIM: %.4f',...
        2^level, psnr_fresh, ssim_fresh))
else
    title('FRESH')
end

% save upsampled result to file

if save_to_file
 imwrite(uint8(img_fresh),[int2str(n) '_4X.png'])
end
end

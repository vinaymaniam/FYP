clear;
I=im2double(imread('92monarch_PSNR35.82.bmp'));
GT=im2double(imread('92monarch.bmp'));
psnr(I(3:end-2,3:end-2),GT(3:end-2,3:end-2))

I=im2double(imread('93zebra584_PSNR32.35.bmp'));
GT=im2double(imread('93zebra584.bmp'));
psnr(I(3:end-2,3:end-2),GT(3:end-2,3:end-2))

I=im2double(imread('94ppt656_PSNR28.77.bmp'));
GT=im2double(imread('94ppt656.bmp'));
psnr(I(3:end-2,3:end-2),GT(3:end-2,3:end-2))

I=im2double(imread('95pepper512_PSNR36.91.bmp'));
GT=im2double(imread('95pepper512.bmp'));
psnr(I(3:end-2,3:end-2),GT(3:end-2,3:end-2))


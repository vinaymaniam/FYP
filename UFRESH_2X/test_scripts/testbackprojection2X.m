clear
dwtmode('spd');
filter = 'db2';
HR = im2double(imread('Testing_Images/GT/Set5/1head280.bmp'));
LR = imresize(HR,0.5);
LR = imresize(LR,2);
SR = LR;
SR = backprojection_2X(SR,HR,filter);
psnr(SR,HR)

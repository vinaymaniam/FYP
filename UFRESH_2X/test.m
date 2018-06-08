clear
y = im2double(imread('Testing_Images/GT/Set5/5childface512.bmp'));

x0 = imresize(imresize(y,0.5),2);
sf = round(size(y)/sqrt(2));
x1 = imresize(imresize(y,sf),size(y));

psnr(x0,y)
psnr(x1,y)
psnr(x0,x1)
psnr(x2,x0)

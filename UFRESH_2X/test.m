clear
y = im2double(imread('Testing_Images/GT/Set5/5childface512.bmp'));

x0 = imresize(y,0.5);
xref = y(1:2:end,1:2:end);

% pp = spline(

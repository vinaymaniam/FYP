function [ img_fresh_shifted ] = fresh_shifted( FRESH, BICU, upscale ) % change the FRESH generated LR interpolated image into standard interpolation LR image
img_bicu=imresize(BICU,(2^upscale),'bicubic');
img_low = imresize(FRESH,1/(2^upscale));
scale = size(FRESH)./size(img_low);
shift = (scale-1)/2;%1.5

[X,Y] = meshgrid((1:size(FRESH,2))+shift(2), (1:size(FRESH,1))+shift(1));
img_fresh_shifted = interp2(FRESH, X, Y);           
img_fresh_shifted(isnan(img_fresh_shifted)) = img_bicu(isnan(img_fresh_shifted)); %Y component

end


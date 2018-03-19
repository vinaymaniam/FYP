function [I_bicshift] = intp_bic(I_low, factor)
    I_bic=imresize(I_low,factor);
    [mm, nn]=size(I_bic); 
    [xx, yy]=meshgrid(1:nn,1:mm); 
    shift=0.5*(factor-1);
    I_bicshift=interp2(xx,yy,I_bic,xx+shift,yy+shift,'spline');
end


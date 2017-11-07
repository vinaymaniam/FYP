function I_bp=backpropogation_2X(I_r,Orig)
% dwtmode('per');
% [Lo_d,Hi_d,Lo_r,Hi_r]=wfilters('bior4.4');
% lf=length(Lo_d);
% sizeEXT = round(length(Lo_d)/2);

I=Orig;
I_up=I_r;

filter='bior4.4';upscale_level=1;
[c_org,l_org] = wavedec2(I,upscale_level,filter);
I_lowc = appcoef2(c_org,l_org,filter,upscale_level)./2^upscale_level;
rangeImg=[min(I_lowc(:)) max(I_lowc(:))];

I_up=range0toN(I_up,rangeImg);
[c1,l1] = wavedec2(I_up,1,filter);
[h1,v1,d1] = detcoef2('all',c1,l1,1);

c_rec=[2*I_lowc(:);  h1(:); v1(:); d1(:);];

[c_org,l_org] = wavedec2(zeros(size(I_up)),1,filter);

I_rec = waverec2(c_rec,l_org,filter);


I_bp=range0toN(I_rec,rangeImg);

end


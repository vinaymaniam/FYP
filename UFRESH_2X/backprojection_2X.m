function I_bp=backprojection_2X(Ir,Orig,filter)  
    [c_org,l_org] = wavedec2(Orig,1,filter);
    Ilowc = appcoef2(c_org,l_org,filter,1)./2;
    
    rangeImg=[min(Ilowc(:)) max(Ilowc(:))];
%     disp(rangeImg)
    Ir=range0toN(Ir,rangeImg);
    [c1,l1] = wavedec2(Ir,1,filter);
    [h1,v1,d1] = detcoef2('all',c1,l1,1);
    
    c_rec=[2*Ilowc(:);  h1(:); v1(:); d1(:);];

    I_rec = waverec2(c_rec,l1,filter);

    I_bp=range0toN(I_rec,rangeImg);
end


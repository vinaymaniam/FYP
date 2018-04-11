function I_bp=backprojection(Ir,Ilowc,filter)        
    rangeImg=[min(Ilowc(:)) max(Ilowc(:))];     
    Ir=range0toN(Ir,rangeImg);    
    [c1,l1] = wavedec2(Ir,1,filter);
    [h1,v1,d1] = detcoef2('all',c1,l1,1);
    
    c_rec=[2*Ilowc(:);  h1(:); v1(:); d1(:);];

    I_rec = waverec2(c_rec,l1,filter);

    I_bp=range0toN(I_rec,rangeImg);
end


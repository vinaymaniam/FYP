function [] = ConvMat2Cell(n)    
    %load(sprintf('pyMap%icell96',n));
    load(sprintf('pyMap%icell96',n));
    Map = reshape(Map,25,25,[]);
    Map = squeeze(num2cell(Map,[1,2]));
    Res = squeeze(num2cell(Res,2));
    save(sprintf('pyMap%icell96',n),'Map','Res');
end
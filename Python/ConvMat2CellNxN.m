function [success] = ConvMat2Cell(n, stage, numneighbours, mode)
    load(sprintf('data_files/pyMap%imat%i',n,numneighbours));
    m = zeros([mode*mode,mode*mode+1,n]);
    for i = 1:n
        m(:,:,i) = Map(i,:,:);
    end
    Map = m;
    Map = squeeze(num2cell(Map,[1,2]));
    Res = squeeze(num2cell(Res,2));
    save(sprintf('data_files/%ipyMap%icell%i_%ix%i',stage,n,numneighbours,mode,mode),'Map','Res');
    success = 1;
end
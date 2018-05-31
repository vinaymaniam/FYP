function [success] = ConvMat2Cell(n, stage, numneighbours)
    load(sprintf('data_files/pyMap%imat%i',n,numneighbours));
    m = zeros([25,26,n]);
    for i = 1:n
        m(:,:,i) = Map(i,:,:);
    end
    Map = m;
    % Map = reshape(Map,25,25,[]);
    Map = squeeze(num2cell(Map,[1,2]));
    Res = squeeze(num2cell(Res,2));
    save(sprintf('data_files/%ipyMap%icell%i',stage,n,numneighbours),'Map','Res');
    success = 1;
end
function [success] = ConvMat2Cell(n, stage, numneighbours)
    load(sprintf('data_files/pyMap%imat%i',n,numneighbours));
    m = zeros([36,37,n]);
    for i = 1:n
        m(:,:,i) = Map(i,:,:);
    end
    Map = m;
    Map = squeeze(num2cell(Map,[1,2]));
    Res = squeeze(num2cell(Res,2));
    save(sprintf('data_files/%ipyMap%icell%i_3x3',stage,n,numneighbours),'Map','Res');
    success = 1;
end
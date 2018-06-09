function [success] = ConvMat2CellL20(n, stage)
    load(sprintf('data_files/pyMap%imat192',n));
    m = zeros([25,26,n]);
    for i = 1:n
        m(:,:,i) = Map(i,:,:);
    end
    Map = m;
    % Map = reshape(Map,25,25,[]);
    Map = squeeze(num2cell(Map,[1,2]));
    Res = squeeze(num2cell(Res,2));
    save(sprintf('data_files/%ipyMap%icell192_L20',stage,n),'Map','Res');
    success = 1;
end
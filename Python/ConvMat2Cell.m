function [success] = ConvMat2Cell(n)
    load(sprintf('data_files/pyMap%imat96',n));
    m = zeros([25,26,n]);
    for i = 1:n
        m(:,:,i) = Map(i,:,:);
    end
    Map = m;
    % Map = reshape(Map,25,25,[]);
    Map = squeeze(num2cell(Map,[1,2]));
    Res = squeeze(num2cell(Res,2));
    save(sprintf('data_files/pyMap%icell96',n),'Map','Res');
    success = 1;
end
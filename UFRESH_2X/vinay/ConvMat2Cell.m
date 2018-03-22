i = 4096;
load(sprintf('pyMap%i',i));
% Map = reshape(Map,25,25,[]);
% Map = squeeze(num2cell(Map,[1,2]));
save(sprintf('pyMap%icell',i),'Map');
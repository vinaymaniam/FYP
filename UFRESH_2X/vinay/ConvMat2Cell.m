i = 8192;
% load(sprintf('pyMap%i',i));
load(sprintf('pyMap%icell96',i));
Map = reshape(Map,25,25,[]);
Map = squeeze(num2cell(Map,[1,2]));
save(sprintf('pyMap%icell96',i),'Map');
i = 4096;
% load(sprintf('pyMap%i',i));
load(sprintf('pyMap%icell96',i));
Map = reshape(Map,25,25,[]);
Map = squeeze(num2cell(Map,[1,2]));
Res = squeeze(num2cell(Res,2));
save(sprintf('pyMap%icell96',i),'Map','Res');
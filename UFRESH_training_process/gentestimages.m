addpath('TrainingData')
addpath('ksvd');
addpath('ksvd/ksvdbox');
addpath('ksvd/ksvdbox/private_ccode');
addpath('ksvd/ompbox');
addpath('utils');
addpath('Step2_Kmeans_clustering');


dirx = 'TrainingData/GT';
pattern = '*.bmp';
xpath = glob(dirx, pattern);

newdir = 'TrainingData/bicubic';
if ~exist(newdir)
    mkdir(newdir);
end
sf = 2;
for i = 1:length(xpath)    
    x = im2double(imread(xpath{i})); 
    y = imresize(imresize(x,0.5),size(x));
    fprintf('%.0f\n',psnr(x,y));
    name = split(xpath{i},'\');
    name = name{end}; 
    imwrite(im2uint8(y),sprintf('%s/%s', newdir, name));
end

% newdir = 'TrainingData/bior44';
% if ~exist(newdir)
%     mkdir(newdir);
% end
% sf = 2;
% filter = 'bior4.4';
% for i = 1:length(xpath)    
%     x = im2double(imread(xpath{i}));
%     [corg,lorg] = wavedec2(x,1,filter);
%     y = appcoef2(corg,lorg,filter,1)./2;  
%     y = imresize(y,size(x));
%     disp([size(x), size(y)])
%     name = split(xpath{i},'\');
%     name = name{end}; 
%     imwrite(im2uint8(y),sprintf('%s/%s', newdir, name));
% end


%% Generate intermediary training images for alternative cascade method

%% STAGE 1.5-2
diry = 'TrainingData/GT';
pattern = '*.bmp';
ypath = glob(diry, pattern);
newdir = 'TrainingData/Alt1';
if ~exist(newdir)
    mkdir(newdir);
end

for i = 1:length(ypath)    
    y = im2double(imread(ypath{i}));
    sf = round(size(y)/sqrt(2));
    x = imresize(imresize(y,sf),size(y));
    x0 = imresize(imresize(y,0.5),2);
    fprintf('%.0f    |     %.0f\n',psnr(x,y),psnr(x0,y))
    name = split(ypath{i},'\');
    name = name{end}; 
    imwrite(im2uint8(x),sprintf('%s/%s', newdir, name));
end

%% STAGE 1-1.5
diry = 'TrainingData/Alt1';
pattern = '*.bmp';
ypath = glob(diry, pattern);
newdir = 'TrainingData/Alt0';
if ~exist(newdir)
    mkdir(newdir);
end

for i = 1:length(ypath)    
    y = im2double(imread(ypath{i}));
    sf = round(size(y)/sqrt(2));
    x = imresize(imresize(y,sf),size(y));
    x0 = imresize(imresize(y,0.5),2);
    fprintf('%.0f    |     %.0f\n',psnr(x,y),psnr(x0,y))
    name = split(ypath{i},'\');
    name = name{end}; 
    imwrite(im2uint8(x),sprintf('%s/%s', newdir, name));
end
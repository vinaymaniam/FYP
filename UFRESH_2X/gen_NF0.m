dirx = 'TrainingData/GT';
pattern = '*.bmp';
ypath = glob(dirx, pattern);
newdir = 'TrainingData/NF_0';
if ~exist(newdir)
    mkdir(newdir);
end
sf = 2;
for i = 1:length(ypath)    
    y = im2double(imread(ypath{i}));
    x = imresize(imresize(y,0.5),2);
    fprintf('%.0f\n',psnr(x,y))
    name = split(ypath{i},'\');
    name = name{end}; 
    imwrite(im2uint8(x),sprintf('%s/%s', newdir, name));
end
dirx = 'TrainingData/GT';
pattern = '*.bmp';
xpath = glob(dirx, pattern);
newdir = 'TrainingData/bior44';
if ~exist(newdir)
    mkdir(newdir);
end
sf = 2;
filter = 'bior4.4';
for i = 1:length(xpath)    
    x = im2double(imread(xpath{i}));
    [corg,lorg] = wavedec2(x,1,filter);
    y = appcoef2(corg,lorg,filter,1)./2;  
    y = imresize(y,size(x));
    disp([size(x), size(y)])
    name = split(xpath{i},'\');
    name = name{end}; 
    imwrite(im2uint8(y),sprintf('%s/%s', newdir, name));
end


% newdir = 'Testing_Images/bicubic/Set5';
% sf = 4;
% for i = 1:length(xpath)
%     x = imread(xpath{i});
%     y = x(1:sf:end,1:sf:end);
%     name = split(xpath{i},'\');
%     name = name{end};    
%     y = imresize(y, sf);
%     disp(psnr(x,y))
%     imwrite(y,sprintf('%s/%s', newdir, name));
% end

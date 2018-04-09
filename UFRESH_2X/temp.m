dirx = 'Testing_Images/GT/Set5';
pattern = '*.bmp';
xpath = glob(dirx, pattern);
newdir = 'Testing_Images/bicubic/Set5';

sf = 4;

for i = 1:length(xpath)
    x = imread(xpath{i});
    y = x(1:sf:end,1:sf:end);
    name = split(xpath{i},'\');
    name = name{end};    
    y = imresize(y, sf);
    disp(psnr(x,y))
    imwrite(y,sprintf('%s/%s', newdir, name));
end

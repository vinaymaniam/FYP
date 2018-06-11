function [ images offsets croppedOriginal ] = SynthDataset(im, numImages, blurSigma)    
    padRatio = 0.2;
    
%     workingRowSub = round(0.5 * padRatio * size(im, 1)) : round((1 - 0.5 * padRatio) * size(im, 1));
%     workingColSub = round(0.5 * padRatio * size(im, 2)) : round((1 - 0.5 * padRatio) * size(im, 2));
    workingRowSub = 1:size(im,1);
    workingColSub = 1:size(im,2);
    croppedOriginal = im(workingRowSub, workingColSub);
%     croppedOriginal = im;

    offsets(1, :) = [ 0 0 ];
    images{1} = im(workingRowSub, workingColSub);    
    
    for i = 2 : numImages
        offsets(i, :) = 2 * rand - 1;
        offsetRowSub = workingRowSub - offsets(i, 2);
        offsetColSub = workingColSub - offsets(i, 1);
        [ x y ] = meshgrid(1 : size(im, 2), 1 : size(im, 1));
        [ x2 y2 ] = meshgrid(offsetColSub, offsetRowSub);
        images{i} = interp2(x, y, im, x2, y2);               
    end
    
    blurKernel = fspecial('gaussian', 3, blurSigma);
    
    for i = 1 : numImages
        images{i} = conv2(images{i}, blurKernel, 'same');
        curIm = images{i};
%         images{i} = curIm(2 : 2 : end - 1, 2 : 2 : end - 1);
        images{i} = curIm(1:2:end, 1:2:end);
    end
end


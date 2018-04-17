function [res] = resLearn(hr, scale, blocksize, heirarchy, index, Map, kernel)
scale = round(log2(scale));

% Artifically downsample the input image
[corg,lorg] = wavedec2(hr, scale, kernel);
lr = appcoef2(corg, lorg, kernel, scale)./2;
lr = imresize(lr, size(hr));
lr = range0toN(lr,[0,1]);

% Super resolve downsampled image
sr = ufresh2(lr, blocksize, heirarchy, index, Map);
sr = range0toN(sr,[0,1]);

% Compute the residue as original - super resolved image
res = hr - sr;
end


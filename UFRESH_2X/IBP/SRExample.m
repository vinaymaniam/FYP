%   Author: Victor May
%   Contact: mayvic(at)gmail(dot)com
%   $Date: 2011-11-19 $
%   $Revision: $
%
% Copyright 2011, Victor May
% 
%                          All Rights Reserved
% 
% All commercial use of this software, whether direct or indirect, is
% strictly prohibited including, without limitation, incorporation into in
% a commercial product, use in a commercial service, or production of other
% artifacts for commercial purposes.     
%
% Permission to use, copy, modify, and distribute this software and its
% documentation for research purposes is hereby granted without fee,
% provided that the above copyright notice appears in all copies and that
% both that copyright notice and this permission notice appear in
% supporting documentation, and that the name of the author 
% not be used in advertising or publicity pertaining to
% distribution of the software without specific, written prior permission.        
%
% For commercial uses contact the author.
% 
% THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO
% THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
% FITNESS FOR ANY PARTICULAR PURPOSE.  IN NO EVENT SHALL THE AUTHOR BE 
% LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL
% DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
% PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
% ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
% THIS SOFTWARE.

% This file is an example script for executing the Super-Resolution
% algorithm.
clc
clear variables
close all
dbstop if error
%% Prepare the reference image
im = imread('rice.png');
% im = rgb2gray(im);
im = imresize(im, [ 128 128 ]);
im = im2double(im);
%% Simulate the low-resolution images
numImages = 4;
blurSigma = 1;
[ images offsets HR ] = SynthDataset(im, numImages, blurSigma);
%% Compute the Super-Resolution image
[ lhs rhs  ] = SREquations(images, offsets, blurSigma);

K = sparse(1 : size(lhs, 2), 1 : size(lhs, 2), sum(lhs, 1));
initialGuess = K \ lhs' * rhs; % This is an 'average' image produced from the LR images.

SR = GradientDescent(lhs, rhs, initialGuess);
SR = reshape(SR, sqrt(numel(SR)), sqrt(numel(SR)));
%% Visualize the results
ShowLRImages(images);
figure;
subplot(1, 3, 1);
imshow(SR);
title('Super-Resolution');
subplot(1, 3, 2);
LR = imresize(imresize(HR, 0.5), size(HR));
imshow(LR);
title('Bicubic Interpolation');
subplot(1, 3, 3);
imshow(HR);
title('Reference');
%% Compute the mean-square error of reconstruction.
mse = sum((SR(:) - HR(:)) .* (SR(:) - HR(:))) / numel(SR);
fprintf(1, 'Reconstruction Mean-square error: %f\n', mse);
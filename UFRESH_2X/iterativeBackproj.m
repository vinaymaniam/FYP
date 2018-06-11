function SR = iterativeBackproj(LR)
%% Simulate the low-resolution images
numImages = 4;
blurSigma = 1;
[ images, offsets , ~ ] = SynthDataset(LR, numImages, blurSigma);
%% Compute the Super-Resolution image
[ lhs, rhs  ] = SREquations(images, offsets, blurSigma);

K = sparse(1 : size(lhs, 2), 1 : size(lhs, 2), sum(lhs, 1));
initialGuess = K \ lhs' * rhs; % This is an 'average' image produced from the LR images.

SR = GradientDescent(lhs, rhs, initialGuess);
SR = reshape(SR, sqrt(numel(SR)), sqrt(numel(SR))); 
end
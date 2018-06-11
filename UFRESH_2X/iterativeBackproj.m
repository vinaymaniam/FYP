%% Only works for square images
function SR = iterativeBackproj(LR)
    % Make input into square image
    sz = min(size(LR));
    orig = LR;
    LR = LR(1:sz,1:sz);
    blurSigma = 1;
    images{1} = LR;
    offsets = [0,0];
    %% Compute the Super-Resolution image
    [ lhs, rhs  ] = SREquations(images, offsets, blurSigma);

    K = sparse(1 : size(lhs, 2), 1 : size(lhs, 2), sum(lhs, 1));
    initialGuess = K \ lhs' * rhs; % This is an 'average' image produced from the LR images.

    SR = GradientDescent(lhs, rhs, initialGuess);
    SR = reshape(SR, size(LR)*2+[1,1]);
    SR = imresize(SR,size(LR)*2);
    % Add in extra bit that was cropped off
    if size(orig,1) > sz
        orig = imresize(orig,2);
        SR = [SR; orig(2*sz+1:end,:)];
    elseif size(orig,2) > sz
        orig = imresize(orig,2);
        SR = [SR, orig(:,2*sz+1:end)];
    end
end
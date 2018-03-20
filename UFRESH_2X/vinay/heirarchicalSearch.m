%% 12s down to 8s by changing xx-2xc+cc to cc-2xc
%% Down to 4s with heirarchical search
%% Down to 3s with cell array in reconstructFromMap.m
%% 12 was for 2048 codewords, 3s is for 4096!!
function [ idx ] = heirarchicalSearch( X,heirarchy )
    idx = zeros(size(X,2),2);
    X = X';
    CC = single(sum(heirarchy(:,:,1).^2, 2)');
    XC = single(X * heirarchy(:,:,1)');
    dists = (bsxfun(@minus, CC, 2*XC));
    [~,idx(:,1)] = min(dists,[],2);
    for i = 1:size(heirarchy,1)
        heir = single(squeeze(heirarchy(i,:,2:end))');
        CC = single(sum(heir.^2, 2)');
        XC = single(X(idx(:,1)==i,:)) * heir';
        dists = (bsxfun(@minus, CC, 2*XC));
        [~,idx2]=min(dists,[],2);
        idx(idx(:,1)==i, 2) = idx2;
    end        
end
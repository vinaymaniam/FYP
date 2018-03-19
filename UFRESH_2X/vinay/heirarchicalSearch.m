%% 12s down to 7s by changing xx-2xc+cc to cc-2xc
function [ idx ] = heirarchicalSearch( X,heirarchy )
    idx = zeros(size(X,2),2);
    X = X';
    CC = single(sum(heirarchy(:,:,1).^2, 2)');
    XC = single(X * heirarchy(:,:,1)');
    dists = (bsxfun(@minus, CC, 2*XC));
    [~,idx(:,1)] = min(dists,[],2);
    
    for i = 1:size(heirarchy,1)
        heir = squeeze(heirarchy(i,:,2:end))';
        CC = single(sum(heir.^2, 2)');
        XC = single(X(idx(:,1)==i,:) * heir');
        dists = (bsxfun(@minus, CC, 2*XC));
        [~,idx2]=min(dists,[],2);
        idx(idx(:,1)==i, 2) = idx2;
    end        
end

%% This function takes the bulk of the time
% Ideally want to cut down here as much as possible. One potential solution
% would be to use a heirarchial search as suggested in 7 Ways
% paper(Timofte)
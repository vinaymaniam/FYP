%% 12s down to 7s by changing xx-2xc+cc to cc-2xc
function [ idx ] = fast_ec( S_all,Center )
S_all=S_all';
Center=Center';

CC = single(sum(Center.^2, 2)');
XC = single(S_all * Center');
% XX = single(sum(S_all.^2, 2));
% dists = (bsxfun(@plus, CC, bsxfun(@minus, XX, 2*XC)));
% Same as dists = CC + XX - 2*XC; but for some reason faster
% Don't need to compute XX for finding minimum
dists = (bsxfun(@minus, CC, 2*XC));

[~,idx]=min(dists,[],2);
end

%% This function takes the bulk of the time
% Ideally want to cut down here as much as possible. One potential solution
% would be to use a heirarchial search as suggested in 7 Ways
% paper(Timofte)
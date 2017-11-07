function [ idx ] = fast_ec( S_all,Center )

S_all=S_all';
Center=Center';

XX = sum(S_all.^2, 2);
CC = sum(Center.^2, 2)';

XC = S_all * Center';
dists = sqrt(bsxfun(@plus, CC, bsxfun(@minus, XX, 2*XC)));
[~,idx]=min(dists,[],2);
end


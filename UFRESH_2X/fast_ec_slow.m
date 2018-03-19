function [ idx ] = fast_ec_slow( S_all,Center )
    idx=[];
    S_all=S_all';
    Center=Center';

    XX = sum(S_all.^2, 2);
    CC = sum(Center.^2, 2)';

    XC = S_all * Center';
    dists = sqrt(bsxfun(@plus, CC, bsxfun(@minus, XX, 2*XC)));
    S=sort(dists,2);
    B=zeros(size(dists));
    for i=1:size(dists,2)
        B(:,i)=S(:,4);
    end
    dists=dists-B;
    for i=1:size(dists,1)
        idx{i}=find(dists(i,:)<0); % idx{i}(1)
    end
end


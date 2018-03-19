function x_rec=linfit(r_est,hi)
N=length(hi);
r_est=[0 r_est N];
x_rec=zeros(N,1);
for ss=1:length(r_est)-1
    idx_org=r_est(ss)+1:r_est(ss+1);
    if length(idx_org)>1
        %% linear fit %%
        idx=idx_org;
        hi_c=hi(:);
        piece=hi_c(idx);
        %     [~,s]=max(abs(piece));
        %     sig=sign(piece(s));
        coefs=polyfit(idx.'-idx(1),piece,1);%%% linear, 1; constance, 
        yfit = coefs(1)*(idx.'-idx(1))+coefs(2);
        x_rec(idx_org)=yfit;
    else
        x_rec(idx_org)=hi(idx_org);
    end
end

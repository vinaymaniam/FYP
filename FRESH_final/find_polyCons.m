function [tt_k]=find_polyCons(K,y_n,lambda,c_m_n_mdf) 

y_n=y_n(:);
N = length(y_n);
T=1/N;
P=N-1;

y_n1 = [y_n(2:end);y_n(end);];
y_d1 = y_n1 - y_n;
s_m=c_m_n_mdf*y_d1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mod(P,2)==1
   J = round((P+1)/2-1);
else
   J = floor((P+1)/2);   
end
uu_k        = acmp_p(s_m, K, J, P, 1);
tt_k        = real(T * log(uu_k) / lambda).';
tt_k        = sort(tt_k);
   

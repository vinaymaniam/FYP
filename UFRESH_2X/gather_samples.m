function [X_cluster,Y_cluster]=gather_samples(Center,idx,X_sample,Y_sample)
clusterszA=64;
X_cluster= [];
Y_cluster= [];
for i=1:3
t=idx(i);
D = pdist2(single(X_sample'),single(Center(:,t)'));
[~, id] = sort(D);                
X_cluster =[X_cluster X_sample(:, id(1:clusterszA))]; 
Y_cluster =[Y_cluster Y_sample(:, id(1:clusterszA))]; 

end


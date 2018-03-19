clear;
load DX_all
load DY_all
load Center1024
cn=size(Center,2);
clusterszA=48; % size of each cluster, can be changed.
Map=[];

% For each centroid Center(:,t)
for t=1:cn
%   Euclidean distance between X and the centroid
%   single just converts to single precision
    D = pdist2(single(X'),single(Center(:,t)'));
%   sort in ascending order, and store the original indices in idx (dv = D(idx))  
    [dv, idx] = sort(D);                
%   Calculate projection matrices and store so they can be called upon at runtime  
    PatchesL = X(:, idx(1:clusterszA)); 
    PatchesH = Y(:, idx(1:clusterszA)); 
    M=PatchesH*PatchesL'*inv(PatchesL*PatchesL'+0.01*eye(size(Center,1))); 
    Map{t}=M;
end
save Map48 Map
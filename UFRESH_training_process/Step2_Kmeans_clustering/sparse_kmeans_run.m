clear;
maxIter =1000;
k  = 1024;
nReplicates    = 1;
load DX_all.mat

tic

[idx,Center,sumd,dist_sketch] = kmeans_sparsified( X, k,'ColumnSamples',true,...
    'Display','iter','Replicates',nReplicates,'Sparsify',false, 'start','++','MaxIter',maxIter);
toc

save Center1024 Center
save idx1024 idx

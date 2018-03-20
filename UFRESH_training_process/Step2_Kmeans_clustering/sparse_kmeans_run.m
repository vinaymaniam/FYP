clear;
maxIter = 300;
nReplicates = 8;
load DX_all.mat

t = []; ks = [];
pool = parpool;
for i = [64 1024 2048 4096 8192]%[256, 1024, 2048, 4096, 8192]
    fprintf('Computing K-means for k=%i...\n',i)
    k = i;
    t1 = tic;        
%     [idx,Center,sumd,dist_sketch] = kmeans_sparsified( X, k,'ColumnSamples',true,...
%         'Display','iter','Replicates',nReplicates,'Sparsify',false, 'start','++','MaxIter',maxIter);
%   Matlab K-means faster    
    stream = RandStream('mlfg6331_64');
    options = statset('UseParallel',1,'UseSubstreams',1,...
    'Streams',stream);
    [idx, Center, ~, ~] = kmeans(X',k,'Options',options,'MaxIter',maxIter,...
                          'Display','final','Replicates',nReplicates);
    
    t = [t, toc(t1)];
    fprintf('K-means for k=%i took %.1f seconds...\n',i,toc(t1))
    ks = [ks, i];
    save(sprintf('Center%i',i), 'Center');
    save(sprintf('idx%i',i), 'idx');
end
figure
plot(ks, t, 'LineWidth', 2)
title('K-means time against # of centroids')
set(gca,'FontSize',20, 'LineWidth',2)
savefig('kmeans_v_time_big_k')

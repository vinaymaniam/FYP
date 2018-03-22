clear;
load DX_all
load DY_all
for i = 4096%[256, 1024, 2048, 4096]
    t1 = tic;
%     load(sprintf('Center%i',i));
    load(sprintf('pyCenter%i',i));
    cn=size(Center,2);
    clusterszA=48; % size of each cluster, can be changed.
    Map=[];
    % For each centroid Center(:,t)
    tsplit = 0;
    for t=1:cn
        t2 = tic;
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
        tsplit = tsplit + toc(t2);
        if mod(t,128) == 0
            fprintf('Finished %i out of %i(%.0f seconds)....\n',t,i,tsplit)
            tsplit = 0;
        end
    end
%     save(sprintf('Map%i',i), 'Map');
    save(sprintf('pyMap%i',i), 'Map');
    fprintf('Mapping Codebook of Size %i Took %.0f seconds....\n\n',i,toc(t1))
    fprintf('-------------------------------------------------\n\n')
end
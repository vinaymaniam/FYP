clear;
load DX_all
load DY_all
for i = 8192%[256, 1024, 2048, 4096]
    t1 = tic;
%     load(sprintf('Center%i',i));
    load(sprintf('pyCenter%i',i));
    cn=size(Center,2);
    clusterszA=96;%96;%48; % size of each cluster, can be changed.
    lambda = 0.01; %Regularisation param, can be changed.
    Map=cell(cn,1);
    % For each centroid Center(:,t)
    for t=1:cn
        D = pdist2(single(X'),single(Center(:,t)'));
    %   sort in ascending order, and store the original indices in idx (dv = D(idx))  
        [dv, idx] = sort(D);                
    %   Calculate projection matrices and store so they can be called upon at runtime  
        PatchesL = X(:, idx(1:clusterszA)); 
        PatchesH = Y(:, idx(1:clusterszA)); 
        M=PatchesH*PatchesL'*inv(PatchesL*PatchesL'+lambda*eye(size(Center,1))); 
        Map{t}=M;   
        fprintf('Finished %i out of %i....\n',t,cn)
    end
%     save(sprintf('Map%i',i), 'Map');
    save(sprintf('pyMap%icell%i',i,clusterszA), 'Map');
    fprintf('Mapping Codebook of Size %i Took %.0f seconds....\n\n',i,toc(t1))
    fprintf('-------------------------------------------------\n\n')
end
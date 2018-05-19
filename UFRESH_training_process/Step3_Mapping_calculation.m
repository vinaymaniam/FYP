clear;
load DX_all
load DY_all
for i = 4096%[256, 1024, 2048, 4096]
    t1 = tic;
%     load(sprintf('Center%i',i));
    load(sprintf('pyCenter%i_',i));
    cn=size(Center,2);
    clusterszA=96;%96;%48; % size of each cluster, can be changed.
    lambdas = [0.01, 0.001, 0.0001];
    tmptime = 0;
    for lambda = 0.01
        Map=cell(cn,1);
        Res = cell(cn,1);
        % For each centroid Center(:,t)
        mypsnr = zeros(cn,1);
        msglen = 0; msg = '';
        for t=1:cn/32
            D = pdist2(single(X'),single(Center(:,t)'));
        %   sort in ascending order, and store the original indices in idx (dv = D(idx)) 
            [dv, idx] = sort(D);   
        %   Calculate projection matrices and store so they can be called upon at runtime  
            a = tic;
            LR = X(:, idx(1:clusterszA)); 
            HR = Y(:, idx(1:clusterszA)); 
            M=HR*LR'*inv(LR*LR'+lambda*eye(size(Center,1))); 
            Map{t}=M;   
            mypsnr(t) = psnr(M*LR,HR);    
            t_iter = toc(t1);
            if mod(t,16) == 0
                msglen = numel(msg);
                msg = sprintf('[Lambda=%.4f] %i | %i [%.0f | %.0fs]',lambda,t,cn,t_iter,time_remaining(t_iter/t,cn-t));  
                fprintf(repmat('\b',1,msglen))
                fprintf(msg)                
            end  
            tmptime = tmptime + toc(a);
        end
        fprintf('\nLambda = %i, PSNR=%.2f\n',lambda,mean(mypsnr))
    end
%     save(sprintf('Map%i',i), 'Map');
    save(sprintf('pyMap%icell%i',i,clusterszA), 'Map');
    fprintf('Mapping Codebook of Size %i Took %.0f seconds....\n\n',i,toc(t1))
    fprintf('-------------------------------------------------\n\n')
end
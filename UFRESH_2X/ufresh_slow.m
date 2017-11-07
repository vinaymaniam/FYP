function [ X_rec_im ] = ufresh_slow( X_test,blocksize,stepsize,Center, Map )
   load DX_1stage; X_sample=single(X); clear X;
   load DY_1stage; Y_sample=single(Y); clear Y;
    X_test_vec = []; % vectorized patches from testing image X;
    cropwidth = size(X_test);
    for j = 1 : stepsize(2) : cropwidth(2)-blocksize(2)+1 % One column is one batch.
		% Rearrange the current batch of image blocks/ patches into corresponding columns block by block . 
		blocks_X = im2colstep(X_test(:, j:j+blocksize(1)-1),blocksize,stepsize);		       
		X_test_vec = [X_test_vec, blocks_X];		
    end
    
    dc_X = mean(X_test_vec);
	X_test_vec = X_test_vec - repmat(dc_X, size(X_test_vec, 1), 1); 
    X_test_vec=single(X_test_vec);
    idx=fast_ec_slow(X_test_vec,Center);
   
    
    X_rec = zeros(size(X_test_vec,1),size(X_test_vec,2)); % recovered X;
    for i=1:size(X_test_vec,2)
        
            [X_cluster,Y_cluster]=gather_samples(Center, idx{i},X_sample,Y_sample);           
            proj=find_project(idx{i},Map,X_cluster, Y_cluster);           
            X_rec(:,i)=proj*X_test_vec(:,i);     
           
    end
    X_rec_mean = X_rec + repmat(dc_X, size(X_rec,1), 1);
    
    %X_rec_mean = X_rec +X_add;
	X_rec_im = col2imstep(X_rec_mean, cropwidth, blocksize, stepsize);
    cnt = countcover(cropwidth,blocksize,stepsize);
	for i = 1:size(cnt,1)
		for j = 1:size(cnt,2)
				if cnt(i,j) == 0
					cnt(i,j) = 1;
				end
		end
    end
    X_rec_im = X_rec_im./cnt;       
end


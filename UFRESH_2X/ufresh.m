function [ X_rec_im ] = ufresh( X_test,blocksize,stepsize,Center, Map )


    X_test_vec = []; % vectorized patches from testing image X;
    cropwidth = size(X_test);
    for j = 1 : stepsize(2) : cropwidth(2)-blocksize(2)+1 % One column is one batch.
		% Rearrange the current batch of image blocks/ patches into 
        % corresponding columns block by block . 
		blocks_X = im2colstep(X_test(:, j:j+blocksize(1)-1),blocksize,stepsize);		       
		X_test_vec = [X_test_vec, blocks_X];		
    end
    
    dc_X = mean(X_test_vec);
	X_test_vec = X_test_vec - repmat(dc_X, size(X_test_vec, 1), 1);
    %% THis line takes the bulk of the time
    idx=fast_ec(X_test_vec,Center);
    %  ------------------------------------
    X_rec = zeros(size(X_test_vec,1),size(X_test_vec,2)); % recovered X;
    for i=1:size(X_test_vec,2)
            s=X_test_vec(:,i);             
            X_rec(:,i)=reshape(Map(:,idx(i)),25,25)*s; 
    end
    X_rec_mean = X_rec + repmat(dc_X, size(X_rec,1), 1);
    
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


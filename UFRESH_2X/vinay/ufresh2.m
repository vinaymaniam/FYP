%% Update from ufresh - Use of heirarchicalKmeans in fast_ec
function [ X_rec_im ] = ufresh2( X_test,blocksize,stepsize,heirarchy,index, Map )


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
    %[heirarchy, index] = heirarchicalKmeans(Center);
    ind=heirarchicalSearch(X_test_vec,heirarchy);
    idx = heir2standard(ind, index);
    %  ------------------------------------
    X_rec_mean = reconstructFromMap(X_test_vec, Map, idx, dc_X);
    
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


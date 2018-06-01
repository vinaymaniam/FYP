%% Update from ufresh2: multiple patch sizes
function [ Xrecim ] = ufresh3( X_test,blocksize,heir5x5,index5x5,Map5x5,heir3x3,index3x3,Map3x3)
    
    cropwidth = size(X_test);
    % vectorized patches from testing image X;
    X_test_vec = im2col(X_test,blocksize,'sliding');    
    %% Split X_test_vec into 2 subproblems: 3x3 and 5x5 based on some metric
    % For now using variance but should probably use something else more
    % representative of high frequency detail
    Xvar = var(X_test_vec);
    % plot(Xvar)
    colidx = Xvar > 0.005;
    Xtest5 = X_test_vec(:,~colidx);
    Xtest3 = X_test_vec(:,colidx);
    %% Need to im2col on Xtest3
    %======================================================================
    dc_X = mean(Xtest5);
	Xtest5 = X_test5 - repmat(dc_X, size(Xtest5, 1), 1);
    %% THis line takes the bulk of the time
    ind=heirarchicalSearch(Xtest5,heir5x5);
    idx = heir2standard(ind, index5x5);
    Xrec5 = reconstructFromMap(X_test_vec, Map5x5, idx, dc_X);
    ind=heirarchicalSearch(Xtest3,heir3x3);
    idx = heir2standard(ind, index3x3);
    Xrec3 = reconstructFromMap(X_test_vec, Map3x3, idx, dc_X);
    Xrecim3 = zeros(size(Xtest5,1),size(Xtest3,2));
    for i = 1:size(Xtest3,2)
        patch5x5 = mergepatch(Xrec3,[3,3],[5,5]);
        Xrecim3(:,i) = im2col(patch5x5,blocksize);
    end
    %  ------------------------------------    
    Xrecmean = zeros(size(X_test_vec));
    Xrecmean(:,~colidx) = Xrec5;
    Xrecim = mergePatch(Xrecmean, blocksize, cropwidth);      
end

function [img] = mergePatch(p, bs, cw)
    y = bs(1); x = bs(2);
    Y = cw(1); X = cw(2);
    img = zeros(Y, X);
    coeff = zeros(Y, X);
    p_idx = 1;
    for xx=1:x
      for yy=1:y
          pp = col2im(p(p_idx,:), [y x], [Y X], 'sliding');
          img(yy:yy+Y-y,xx:xx+X-x) = img(yy:yy+Y-y,xx:xx+X-x)+pp;
          coeff(yy:yy+Y-y,xx:xx+X-x) = coeff(yy:yy+Y-y,xx:xx+X-x)+1;
          p_idx = p_idx+1;
      end
    end
    img = img ./ coeff;
end




% OLD
% function [ X_rec_im ] = ufresh2( X_test,blocksize,heirarchy,index, Map )
% 
%     X_test_vec = []; % vectorized patches from testing image X;
%     cropwidth = size(X_test);
%     for j = 1:cropwidth(2)-blocksize(2)+1 % One column is one batch.
% 		% Rearrange the current batch of image blocks/ patches into corresponding columns block by block . 
% 		blocks_X = im2colstep(X_test(:, j:j+blocksize(1)-1),blocksize,[1,1]);		       
% 		X_test_vec = [X_test_vec, blocks_X];		
%     end
%     
%     dc_X = mean(X_test_vec);
% 	X_test_vec = X_test_vec - repmat(dc_X, size(X_test_vec, 1), 1); 
% 
%     ind=heirarchicalSearch(X_test_vec,heirarchy);
%     idx = heir2standard(ind, index);
%     %  ------------------------------------
%     X_rec_mean = reconstructFromMap(X_test_vec, Map, idx, dc_X);
%     
% 	X_rec_im = col2imstep(X_rec_mean, cropwidth, blocksize, [1,1]);
%     cnt = countcover(cropwidth,blocksize,[1,1]);
% 	for i = 1:size(cnt,1)
% 		for j = 1:size(cnt,2)
% 				if cnt(i,j) == 0
% 					cnt(i,j) = 1;
% 				end
% 		end
%     end
%     X_rec_im = X_rec_im./cnt;       
% end


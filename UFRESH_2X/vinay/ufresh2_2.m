%% Update from ufresh - Use of heirarchicalKmeans in fast_ec
function [ Xrecim ] = ufresh2_2( X_test,blocksize,heirarchy,index, Map )
    
    cropwidth = size(X_test);
    % vectorized patches from testing image X;
    X_test_vec = im2col(X_test,blocksize,'sliding');
    dc_X = mean(X_test_vec);
	X_test_vec = X_test_vec - repmat(dc_X, size(X_test_vec, 1), 1);
    %% THis line takes the bulk of the time
    idx=heirarchicalSearch_2(X_test_vec,heirarchy, index);
%     idx = heir2standard(ind, index);
    %  ------------------------------------
    Xrecmean = reconstructFromMap(X_test_vec, Map, idx, dc_X);
    
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


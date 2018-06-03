%% Update from ufresh2: multiple patch sizes
% 6x6 and 3x3 <-- 3x3 non overlapping
function [ Xrecim ] = ufresh4( X_test,blocksize,heir6x6,index6x6,Map6x6,heir3x3,index3x3,Map3x3)
    
    cropwidth = size(X_test);
    % vectorized patches from testing image X;
    X_test_vec = im2col(X_test,blocksize,'sliding');    
    %% Split X_test_vec into 2 subproblems: 3x3 and 5x5 based on some metric
    % For now using variance but should probably use something else more
    % representative of high frequency detail
    Xvar = var(X_test_vec);
    % plot(Xvar)
    colidx = Xvar > 0.01;
    Xtest6 = X_test_vec(:,~colidx);
    dc_X6 = mean(Xtest6);
	Xtest6 = Xtest6 - repmat(dc_X6, size(Xtest6, 1), 1);   
    %% Need to im2col on Xtest3
    Xt3 = X_test_vec(:,colidx);
    if size(Xt3,2) > 0
        a  = reshape(Xt3,6,6,[]);
        Xtest3 = [];
        for i = 1:size(a,3)
            Xtest3(:,1+4*(i-1):4+4*(i-1)) = im2col(a(:,:,i),[3,3],'distinct');
        end
        dc_X3 = mean(Xtest3);
        Xtest3 = Xtest3 - repmat(dc_X3, size(Xtest3, 1), 1);
    end
    %======================================================================    
    %% THis line takes the bulk of the time
    ind=heirarchicalSearch(Xtest6,heir6x6);
    idx = heir2standard(ind, index6x6);
    Xrec6 = reconstructFromMap(Xtest6, Map6x6, idx, dc_X6);
    if size(Xt3,2) > 0
        ind=heirarchicalSearch(Xtest3,heir3x3);
        idx = heir2standard(ind, index3x3);
        Xrec3 = reconstructFromMap(Xtest3, Map3x3, idx, dc_X3);
        Xrecim3 = zeros(size(Xtest6,1),size(Xtest3,2)/4);
        for i = 1:(size(Xtest3,2)/4)
            patch6x6 = mergePatch(Xrec3(:,1+9*(i-1):9+9*(i-1)),[3,3],[5,5]);
            patch6x6 = col2im(Xrec3(:,1+4*(i-1):4+4*(i-1)),[3,3],[6,6]);
            Xrecim3(:,i) = reshape(patch6x6,36,1);
        end
    end
    %  ------------------------------------    
    Xrecmean = zeros(size(X_test_vec));
    Xrecmean(:,~colidx) = Xrec6;
    Xrecmean(:,colidx) = Xrecim3;
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



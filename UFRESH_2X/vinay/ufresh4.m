%% Update from ufresh2: multiple patch sizes
% 6x6 and 3x3 <-- 3x3 non overlapping
function [ Xrecim ] = ufresh4( X_test,blocksize,heirN1,indexN1,MapN1,heirN2,indexN2,MapN2)
    
    cropwidth = size(X_test);
    % vectorized patches from testing image X;
    X_test_vec = im2col(X_test,blocksize,'sliding');    
    %% Split X_test_vec into 2 subproblems: N1xN1 and N2xN2 based on some metric
    % For now using variance but should probably use something else more
    % representative of high frequency detail
    Xvar = var(X_test_vec);
    % plot(Xvar)
    colidx = Xvar > 0.01;
    XtestN1 = X_test_vec(:,~colidx);
    dc_XN1 = mean(XtestN1);
	XtestN1 = XtestN1 - repmat(dc_XN1, size(XtestN1, 1), 1);   
    %% Need to im2col on Xtest3
    XtN2 = X_test_vec(:,colidx);
    if size(XtN2,2) > 0
        a  = reshape(XtN2,6,6,[]);
        XtestN2 = [];
        for i = 1:size(a,3)
            XtestN2(:,1+4*(i-1):4+4*(i-1)) = im2col(a(:,:,i),[3,3],'distinct');
        end
        dc_X3 = mean(XtestN2);
        XtestN2 = XtestN2 - repmat(dc_X3, size(XtestN2, 1), 1);
    end
    %======================================================================    
    %% THis line takes the bulk of the time
    ind=heirarchicalSearch(XtestN1,heirN1);
    idx = heir2standard(ind, indexN1);
    XrecN1 = reconstructFromMap(XtestN1, MapN1, idx, dc_XN1);
    if size(XtN2,2) > 0
        ind=heirarchicalSearch(XtestN2,heirN2);
        idx = heir2standard(ind, indexN2);
        XrecN2 = reconstructFromMap(XtestN2, MapN2, idx, dc_X3);
        XrecimN2 = zeros(size(XtestN1,1),size(XtestN2,2)/4);
        for i = 1:(size(XtestN2,2)/4)           
            patch6x6 = col2im(XrecN2(:,1+4*(i-1):4+4*(i-1)),[3,3],[6,6],'distinct');
            XrecimN2(:,i) = reshape(patch6x6,36,1);
        end
    end
    %  ------------------------------------    
    Xrecmean = zeros(size(X_test_vec));
    Xrecmean(:,~colidx) = XrecN1;
    Xrecmean(:,colidx) = XrecimN2;
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



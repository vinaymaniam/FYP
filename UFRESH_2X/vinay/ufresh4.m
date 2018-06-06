%% Update from ufresh2: multiple patch sizes
% 6x6 and 3x3 <-- 3x3 non overlapping
function [ Xrecim ] = ufresh4( X_test,blocksize,heirN1,indexN1,MapN1,heirN2,indexN2,MapN2)
    N1 = sqrt(size(heirN1,2));
    N2 = sqrt(size(heirN2,2));
    cropwidth = size(X_test);
    % vectorized patches from testing image X;
    X_test_vec = im2col(X_test,blocksize,'sliding');    
    %% Split X_test_vec into 2 subproblems: N1xN1 and N2xN2 based on some metric
    % For now using variance but should probably use something else more
    % representative of high frequency detail
    Xvar = var(X_test_vec);
    % plot(Xvar)
    colidx = Xvar > 0.025;
    XtestN1 = X_test_vec(:,~colidx);
    dc_XN1 = mean(XtestN1);
	XtestN1 = XtestN1 - repmat(dc_XN1, size(XtestN1, 1), 1);   
    %% Need to im2col on Xtest3
    XtN2 = X_test_vec(:,colidx);
    if size(XtN2,2) > 0
        a = reshape(XtN2, N1, []);
        XtestN2 = im2col(a,[N2,N2],'distinct');
        dc_XN2 = mean(XtestN2);
        XtestN2 = XtestN2 - repmat(dc_XN2, size(XtestN2, 1), 1);
    end
    %======================================================================    
    Xrecmean = zeros(size(X_test_vec));
    %% THis line takes the bulk of the time
    ind=heirarchicalSearch(XtestN1,heirN1);
    idx = heir2standard(ind, indexN1);
    XrecN1 = reconstructFromMap(XtestN1, MapN1, idx, dc_XN1);
    if size(XtN2,2) > 0
        ind=heirarchicalSearch(XtestN2,heirN2);
        idx = heir2standard(ind, indexN2);
        XrecN2 = reconstructFromMap(XtestN2, MapN2, idx, dc_XN2);
        tmp = col2im(XrecN2,[N2,N2],[N1,N1*size(XtestN2,2)/4],'distinct');
        XrecimN2 = reshape(tmp,N1*N1,[]);
        Xrecmean(:,colidx) = XrecimN2;
    end
    %  ------------------------------------        
    Xrecmean(:,~colidx) = XrecN1;    
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



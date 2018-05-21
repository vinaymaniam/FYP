%% Added bias term
function Xrecmean = reconstructFromMap(X_test_vec, Map, idx, dc_X)
    Xrec = zeros(size(X_test_vec,1),size(X_test_vec,2)); % recovered X;
    for i=1:size(X_test_vec,2)
        s=X_test_vec(:,i);
        Xrec(:,i)=Map{idx(i)}*[s; ones(1,size(s,2))];
    end
    Xrecmean = Xrec + repmat(dc_X, size(Xrec,1), 1);
end
%% Using map as cell array MUCH faster than using multidim array
% function Xrecmean = reconstructFromMap(X_test_vec, Map, idx, dc_X)
%     Xrec = zeros(size(X_test_vec,1),size(X_test_vec,2)); % recovered X;
%     for i=1:size(X_test_vec,2)
%         s=X_test_vec(:,i);
%         Xrec(:,i)=Map{idx(i)}*s;
%     end
%     Xrecmean = Xrec + repmat(dc_X, size(Xrec,1), 1);
% end
%% Original version
% function Xrecmean = reconstructFromMap(X_test_vec, Map, idx, dc_X)
%     Xrec = zeros(size(X_test_vec,1),size(X_test_vec,2)); % recovered X;
%     for i=1:size(X_test_vec,2)
%             s=X_test_vec(:,i);             
%             Xrec(:,i)=reshape(Map(:,idx(i)),25,25)*s; 
%     end
%     Xrecmean = Xrec + repmat(dc_X, size(Xrec,1), 1);
% end
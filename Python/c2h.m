% function [heirarchy, index] = c2h(center)
%     % Do center to heirarchy but instead of choosing num neighbours we just
%     % assign by label
%     k = ceil(sqrt(size(center,2)));
%     [labels, h] = kmeans(center',k,'MaxIter',500,'Replicates',3);
%     h = h';
%     for i = 1:size(h,2)
%         heirarchy{i} = center(:, labels==i);
%         index{i} = find(labels==i);
%     end
%     heirarchy{size(h,2)+1} = single(h);
% end
% 
function [heirarchy, index] = c2h(center)
    % Do center to heirarchy but instead of choosing num neighbours we just
    % assign by label
    k = ceil(sqrt(size(center,2)));
    [labels, h] = kmeans(center',k,'MaxIter',500,'Replicates',3);
    h = h';
    for i = 1:size(h,2)
        idx = find(labels==i);
%         knn = knnsearch(center',h(:,i)','K',2*length(idx));
        heirarchy{i} = center(:, labels==i);
        index{i} = find(labels==i);
    end
    heirarchy{size(h,2)+1} = single(h);
end


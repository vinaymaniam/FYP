function idx = heir2standard(ind, index)
% Convert heirarchical indices back to regular indices so they can be used
% to access regressions in 'Map'
    for i = 1:size(ind,1)
        idx(i) = index(ind(i,1),ind(i,2));
    end
end
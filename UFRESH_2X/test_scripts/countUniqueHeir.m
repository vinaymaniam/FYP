%% Checks how many centroids are actually in the heirarchy to see how many
%  slipped through the cracks
clear;
stage = 1;


for n = [128,256,512,1024,2048,4096,8192,16384]
    load(sprintf('%ipyHeirarchy%i_NF',stage,n));
    h = single(heirarchy);   
    H = zeros(size(h,1)*(size(h,3)-1),size(h,2));
    for i = 0:size(h,1)-1
        H((size(h,3)-1)*i+1:(size(h,3)-1)*(i+1),:) = squeeze(h(i+1,:,2:end))';
    end

    a = unique(H,'rows');
    fprintf('%.1f %% missing centroids\n',100*(n-size(a,1))/n)    
    missing(log2(n)-6) = 100*(n-size(a,1))/n;
end
plot(missing)

% CONCLUSION- not all centroids are actually stored in here
% SOLUTION- instead of storing sam number of centroids for each
% hierarchical centroid, just use the labels so some will have more
% children than others
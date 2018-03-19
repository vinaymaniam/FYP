function [heirarchy, ind] = heirarchicalKmeans(centroids)
% centroids is the codebook of existing centroids from standard kmeans
% n is the number of child centroids to associate with each parent centroid
centroids = centroids';
nctrds = round(sqrt(size(centroids,1)));
[~, C] = kmeans(centroids, nctrds);

idx = knnsearch(centroids, C, 'K', 3*nctrds, 'Distance', 'euclidean');

heirarchy = zeros(nctrds, 25, size(idx,2)+1);
heirarchy(:,:,1) = C;
ind = zeros(nctrds, size(idx,2));
for i = 1:size(idx,2)
    for j = 1:nctrds
        heirarchy(j,:,i+1) = centroids(idx(j,i),:);
        ind(j,i) = idx(j,i);
    end
end
end
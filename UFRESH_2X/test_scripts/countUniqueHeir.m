%% Checks how many centroids are actually in the heirarchy to see how many
%  slipped through the cracks
clear;
stage = 1;

uv = [];
ns = [128,256,512,1024,2048,4096,8192,16384];
for n = ns
    load(sprintf('%ipyHeirarchy%i_NF',stage,n));    
    for k = 1:4
        endpoint = int32(k*(size(heirarchy,3)-1)/4);
        h = single(heirarchy(:,:,1:endpoint));   
        H = zeros(size(h,1)*(size(h,3)-1),size(h,2));       
        for i = 0:size(h,1)-1
            H((size(h,3)-1)*i+1:(size(h,3)-1)*(i+1),:) = squeeze(h(i+1,:,2:end))';
        end
        a = unique(H,'rows');
        uv(log2(n)-6,k) = 100*(n-size(a,1))/n;        
    end    
end
missing = mean(uv,1);
figure;
plotme(1:length(missing),uv,'mix',1)
title('Percentage of Centroids Lost as a Function of \alpha')
xlabel('\alpha')
ylabel('% Lost')
axis([-inf,inf,0,inf])
legend(cellstr(num2str(ns')))
set_gca

% CONCLUSION- not all centroids are actually stored in here
% SOLUTION- instead of storing sam number of centroids for each
% hierarchical centroid, just use the labels so some will have more
% children than others
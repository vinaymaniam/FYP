function proj=find_project(idx,Map,X_cluster, Y_cluster)
Map=single(Map);
error=zeros(3,1);
P=[];
for i=1:3
P{i}=reshape(Map(:,idx(i)),25,25);
res=P{i}*X_cluster-Y_cluster;
error(i)=sum(sum(res.^2));
end
rn=find(error==min(error));
proj=P{rn};
end


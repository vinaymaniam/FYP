function [ I_rec ] = reconst(I)
load Center2048
load Map2048
p=5;
[HI,WI]=size(I); 
C=zeros(p*p,(HI-p+1),(WI-p+1));
I_rec=zeros(size(I));
count=I_rec;
for i=1:1:HI-p+1
    for j=1:1:WI-p+1
    C(:,i,j)=im2col(I(i:i+p-1,j:j+p-1),[p p]);
    end
end

for M=1:HI-p+1   
    for N=1:WI-p+1
        C_patch=C(:,M,N)-mean(C(:,M,N));
        if sum(C_patch.^2, 1)>0.1 % threshold
        % find the projection
        MSE_rough=sqrt(sum((Center-repmat(C_patch,1,size(Center,2))).^2));
        mse=sort(MSE_rough(:));
        t=find(MSE_rough<=mse(1),1); % find the most minimum MSE       
        avgp=Map{t};       
        % multiply by the projection                    
        patch_rec=avgp*C_patch+mean(C(:,M,N));
        else
        patch_rec=C(:,M,N);
        end
        patch_rec=reshape(patch_rec,[p p]);
        % assemble into the reconstructed image
        I_rec(M:M+p-1,N:N+p-1)=I_rec(M:M+p-1,N:N+p-1)+patch_rec;
        count(M:M+p-1,N:N+p-1)= count(M:M+p-1,N:N+p-1)+1; 
    end
end
I_rec=I_rec./count;
end



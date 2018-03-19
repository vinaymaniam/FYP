function [I_rec, IDX, MSE]=FRI_errRepair_ann_lm(I,I_fri,I_fri2,level)
size_I=size(I_fri2);
IDX=[];MSE=[];
ratio=size(I_fri2,2)./size(I_fri,2);
 if level>=2
	patchsize=[5 5];
	num=12; % default 8
 else
    patchsize=[5 5];
    num=12;%default 10
end
p_size=floor(patchsize(1)./2);
I = padarray(I,[num num]);
I_fri = padarray(I_fri,[num num]);
I_fri2 = padarray(I_fri2,[ceil(num*ratio) ceil(num*ratio)]);
range=num-1;
patch_count=0;
I_rec=zeros(size(I_fri2));count=I_rec;
for xx=ceil(num*ratio)+p_size+1:size_I(1)+ceil(num*ratio)-p_size
    fprintf('%i out of %i\n', xx, size_I(1)+ceil(num*ratio)-p_size);
    for yy=ceil(num*ratio)+p_size+1:size_I(2)+ceil(num*ratio)-p_size
        patch_count=patch_count+1;
        patch_fri2_current=I_fri2(xx+(-p_size:p_size),yy+(-p_size:p_size));
        patch_fri2_current=patch_fri2_current(:);
        Region_fri2=I_fri2(round(xx-range*ratio):round(xx+range*ratio),round(yy-ratio*range):round(yy+ratio*range));
        Region_fri=I_fri(round(xx/ratio-range):round(xx/ratio+range),round(yy/ratio-range):round(yy/ratio+range));
        Region_I=I(round(xx/ratio-range):round(xx/ratio+range),round(yy/ratio-range):round(yy/ratio+range));

        Patch_fri=im2col(Region_fri,patchsize,'sliding');

        patch_MSE=sqrt(sum((Patch_fri-repmat(patch_fri2_current,1,size(Patch_fri,2))).^2));
        [mse,idx]=sort(patch_MSE);
        mse=mse(1:4);idx=idx(1:4);

        IDX{patch_count}= idx(:);
        MSE{patch_count}= mse(:);

        [x_p,y_p]=ind2sub(size(Region_fri)-patchsize(1)+1,idx);
        Y_l=[];Y_h=[];
        for ss=1:length(x_p)
            patch_fri=Region_fri(x_p(ss)+(0:patchsize(1)-1),y_p(ss)+(0:patchsize(2)-1));
            patch_I=Region_I(x_p(ss)+(0:patchsize(1)-1),y_p(ss)+(0:patchsize(2)-1));
            Y_l=[Y_l patch_fri(:)];
            Y_h=[Y_h patch_I(:)];  

        end

        YY=Y_l*Y_l';
        M=Y_h*Y_l'*inv(YY+0.1*eye(size(YY)));
        patch_up=M*patch_fri2_current(:);
        patch_up=reshape(patch_up,patchsize);
        I_rec(xx+(-p_size:p_size),yy+(-p_size:p_size))=I_rec(xx+(-p_size:p_size),yy+(-p_size:p_size))+patch_up;%I_fri2(xx+(0:patchsize(1)-1),yy+(0:patchsize(2)-1))+patch_err;
        count(xx+(-p_size:p_size),yy+(-p_size:p_size))=count(xx+(-p_size:p_size),yy+(-p_size:p_size))+1;
    end
end
count(count<1)=1;
I_rec=I_rec./count;
I_rec=I_rec(ceil(num*ratio)+1:end-ceil(num*ratio),ceil(num*ratio)+1:end-ceil(num*ratio));
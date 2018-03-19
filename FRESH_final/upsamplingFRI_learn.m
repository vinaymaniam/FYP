% -------------------------------------------------------------------------
% Communications and Signal Processing Group
% Department of Electrical and Electronic Engineering
% Imperial College London, 2015
%
% Date        : 12/09/2015
% Supervisor  : Dr Pier Luigi Dragotti
% Author      : Xiaoyao Wei
% -------------------------------------------------------------------------
name_vec={'butterfly256'};%{'comic360','butterfly256','bird288','woman344','pepper512','lena512','','childface512','zebra584'};
filter_vec={'bior4.4'};%{'rbio2.8'};
for image_count=1:length(name_vec)
    clearvars -except image_count name_vec filter_vec;
    close all
    tic
    %%%%%%%%%%%%%%%
    atf=1;%1 means artificial downsampling and then perform upsampling. 0 means perform upsampling of photos we took
    %%%%%%%%%%%%%%
    b=1.6;%intermediate level, scale of 1.6
    upscale_level=2; %%% upsampling factor: 2^upscale_level
    lm_iter=0;upres=8;
    diag_iter=3;
    if atf==1
        filter=filter_vec{image_count};filter_fri=filter;a=1;
    else
        filter='rbio6.8';filter_fri='spln6';a=1.3;
    end
    dwtmode('per');
    [Lo_d,Hi_d,Lo_r,Hi_r] = wfilters(filter);
    lf=length(Lo_d);
    sizeEXT = round(length(Lo_d)/2);


    if atf==1
        % % input high-resolution image%% 
        imagename=name_vec{image_count};
        I=double((imread([imagename '.bmp'])));
%         I=double((imread([imagename '.jpg'])));
        [c_org,l_org] = wavedec2(I,upscale_level,filter_fri);
        I_lowc = appcoef2(c_org,l_org,filter,upscale_level)./2^upscale_level;
        %%%%%%%%
    else
        I_lowc=(imread('gate_400d.JPG')); %gate_400d, roadsign_400d
        %I_lowc=I_lowc(1120+(1:100),1675+(1:100),:);%va
        %I_lowc=I_lowc(1130+(1:100),1220+(1:100),:);%pillar
        %I_lowc=I_lowc(1050+50+(1:64),1590+((1:64)),:);%roadsign
        I_lowc=I_lowc(1115+(1:64),1663+(1:64),:);%gate
        %I_lowc=imread('child_small.png');
    end

    [N1,N2,~]=size(I_lowc);NN1=2^upscale_level*N1;NN2=2^upscale_level*N2;
    I_zero=zeros(NN1,NN2);
    [c_org,l_org] = wavedec2(I_zero,upscale_level,filter_fri);
    layer=size(I_lowc,3);
    if layer==1
        I_lowy=I_lowc;
    else
        I_lowy = rgb2ycbcr(I_lowc);
    end
    rangeImg=[min(I_lowy(:)) max(I_lowy(:))];
    I_lin_c=(zeros(N1*2^upscale_level,N2*2^upscale_level,layer));
    for ss=1:layer
        c=zeros(size(c_org));
        cA=I_lowy(:,:,ss);
        c(1:N1*N2)=cA(:);
        c(N1*N2+1:end)=0;
        I_lin_c(:,:,ss)= intp_bic(double(cA),2^upscale_level);
    end
    I_low=double(I_lowy(:   ,:,1)).*2^upscale_level;

    I_int=I_low./2^upscale_level;
    I_org{1}=I_int;


    I_lin = I_lin_c(:,:,1).*2^upscale_level;
    if layer == 3
        hImcb = I_lin_c(:,:,2);
        hImcr = I_lin_c(:,:,3);
    end
    I_lin_org=I_lin;
    I_f=I_lin;
    I_bic=intp_bic(I_int,2^upscale_level);
    up_flag=1; 
    level_up=0;
    level=1;

    for it=1:upscale_level
        if it>1
           I_org{it}=double(I_upres{it-1});
        end
        for supres=0:2+lm_iter
            level=1;
            if supres==0 %downsampling the current high-resolution version and upsampling by FRI
                I_temp=intp_bic(I_org{it},b);
                [c_org,l_org] = wavedec2(I_temp,level,filter_fri);
                I_low = appcoef2(c_org,l_org,filter_fri,level); % the low resolution image we have access to
                [N1,N2]=size(I_low);[NN1,NN2]=size(I_temp);
                c_org(N1*N2+1:end)=0;
                I_lin = waverec2(c_org,l_org,filter); % linear reconstruction of the low resolution image
                I_f=I_lin;
            elseif supres==1 % upsample the current high-resolution version
                I_low=I_org{it}*2;
                [N1,N2]=size(I_low);NN1=2^level*N1;NN2=2^level*N2;
                [c_org,l_org] = wavedec2(zeros(NN1,NN2),level,filter_fri);
                c_org(1:N1*N2)=I_low(:);
                c_org(N1*N2+1:end)=0;
                I_lin = waverec2(c_org,l_org,filter); % linear reconstruction of the low resolution image
                I_lin=range0toN(I_lin,rangeImg);
                I_f=I_lin;
            else
                [c_org,l_org] = wavedec2(I_up1to2{it},level,filter_fri);
                I_low = appcoef2(c_org,l_org,filter_fri,level); % the low resolution image we have access to
                [N1,N2]=size(I_low);[NN1,NN2]=size(I_up1to2{it});
                c_org(N1*N2+1:end)=0;
                I_lin = waverec2(c_org,l_org,filter); % linear reconstruction of the low resolution image
                I_f=I_lin;
            end
            %%%%%%%%%  linear reconstruction along columns, will be used for
            %%%%%%%%%  calculating FRI image Coef_x

            intp_x = I_low;
            for ss=1:level
                intp_x = dyadup(intp_x,'row',0,1);
                intp_x = wextend('addrow','per',intp_x,sizeEXT);
                intp_x = conv2(intp_x' ,Lo_r(:)','full');
                intp_x = intp_x';
                intp_x = intp_x(lf:lf+N1*(2^(ss))-1,:);
            end

            %%%%%%%%%  linear reconstruction along rows, will be used for
            %%%%%%%%%  calculating FRI image Coef_y
            intp_y = I_low;
            for ss=1:level
                intp_y = dyadup(intp_y,'col',0,1);
                intp_y = wextend('addcol','per',intp_y,sizeEXT);
                intp_y = conv2(intp_y ,Lo_r(:)','full');
                intp_y = intp_y(:,lf:lf+N2*(2^(ss))-1);
            end

            kk_ini=floor(size(I_int,1)./4);
            if kk_ini>25
               kk_ini=25;
            end
            kk=kk_ini+mod(supres,2)*8+(it-1)*15;%% no. of discontinuities in one line 

            [phi,~,~,~,xval] = wavefun(filter_fri,level);
            tmp1=find(phi~=0,1);
            tmp2=find(phi~=0,1,'last');
            phi=phi(tmp1:tmp2);
            loc_shift=floor((length(phi)-1)/2);

            %[Corg,Sorg]=wavedec(zeros(NN2,1),level,filter_fri);
            [C_up,S_up]=wavedec(zeros(NN2*(upres),1),log2(upres)+1,filter_fri);
            [t_org,~,~]=setup_filter(filter_fri,N2,level,a);  
            [t,lambda,c_m_n_mdf]=setup_filter(filter_fri,N2,level,a);  T_s=t(2)-t(1);
            %%%%%% upsampling line by line using FRI : horizontal line %%%%%%%%
            Coef_x=[];
            for ss=1:NN1

                y_n=intp_x(ss,1:N2);
                x=zeros(1,NN2);
                x_f=I_f(ss,:);

                inloc_est=find_polyCons(kk,y_n,lambda,c_m_n_mdf)-loc_shift*T_s;
                r_upres=round((inloc_est-t_org(1))/(t_org(2)-t_org(1))*upres+1);
                r_upres=sort(r_upres((r_upres>1)&(r_upres<NN2*upres)));
                rr_est=unique(r_upres);
                lin_rec_up=waverec([y_n(:).' zeros(1,length(C_up)-S_up(1))],S_up,filter_fri);   
                x_est_up=linfit(rr_est,lin_rec_up);
                x_rec_up=wavedec(x_est_up,log2(upres),filter_fri);
                x_rec=x_rec_up(1:NN2);
                Coef_x=[Coef_x;x_rec(:).'];
            end

            [Corg,Sorg]=wavedec(zeros(NN1,1),level,filter_fri);
            [C_up,S_up]=wavedec(zeros(NN1*(upres),1),log2(upres)+1,filter_fri);
            [t_org,~,~]=setup_filter(filter_fri,N1,level,a);  
            [t,lambda,c_m_n_mdf]=setup_filter(filter_fri,N1,level,a);T_s=t(2)-t(1);
            Coef_x=range0toN(Coef_x,rangeImg);
            %%%%%%% upsampling line by line using FRI : vertical line %%%%%%%%
            Coef_y=[];
            for ss=1:NN2
                y_n=intp_y(:,ss);
                x=zeros(NN1,1);
                x_f=I_f(:,ss);

                inloc_est=find_polyCons(kk,y_n,lambda,c_m_n_mdf)-loc_shift*T_s;
                    r_upres=round((inloc_est-t_org(1))/(t_org(2)-t_org(1))*upres+1);
                    r_upres=sort(r_upres((r_upres>1)&(r_upres<NN1*upres)));
                    rr_est=unique(r_upres);
                    lin_rec_up=waverec([y_n(:).' zeros(1,length(C_up)-S_up(1))],S_up,filter_fri);   
                    x_est_up=linfit(rr_est,lin_rec_up);
                    x_rec_up=wavedec(x_est_up,log2(upres),filter_fri);
                    x_rec=x_rec_up(1:NN1);
                Coef_y=[Coef_y x_rec(:)];
            end

            Coef_y=range0toN(Coef_y,rangeImg);


            % %%% add the high freq coefficients estimated using FRI method to the linear reconstruction
            % 
            %dwtmode('per');
            if (supres==1)
                level=it;
            end

            [c1,l1] = wavedec2(Coef_y,level,filter);
            [c2,l2] = wavedec2(Coef_x,level,filter);  
            [c3,l3] = wavedec2(Coef_y,level,filter);

            [h1,~,~] = detcoef2('all',c1,l1,1);
            [~,v1,~] = detcoef2('all',c2,l2,1);
            [~,~,d1] = detcoef2('all',c3,l3,1);
            [c_org,l_org] = wavedec2(zeros(NN1,NN2),level,filter);
            if (it==2)&&(supres==1)
                [h2,~,~] = detcoef2('all',c1,l1,2);
                [~,v2,~] = detcoef2('all',c2,l2,2);
                [~,~,d2] = detcoef2('all',c3,l3,2);
                c_rec=[I_int(:)*4; h2(:); v2(:); d2(:); h1(:); v1(:); d1(:);];
            else
                c_rec=[I_low(:); h1(:); v1(:); d1(:);];
            end
            I_rec = waverec2(c_rec,l_org,filter);
            I_rec=range0toN(I_rec,rangeImg);

            I_f=I_rec;
            I_xy=I_f;

            for iter=1:diag_iter %%%%%% diagnal upsampling (may be not necessary)
                level=1;
                [phi,~,~,~,xval] = wavefun(filter_fri,level);
                % figure,plot(xval,phi)
                tmp1=find(phi~=0,1);
                tmp2=find(phi~=0,1,'last');
                phi=phi(tmp1:tmp2);
                loc_shift=ceil((length(phi)-1)/2);

                  %%%%%%% first downsampling and then upsampling line by line using FRI : diagonal line 45 degree %%%%%%%%
                I_d=I_f;
                Coef_d1=zeros(NN1,NN2);N1=round(NN1./2^level);N2=round(NN2./2^level);
                count=zeros(NN1,NN2);
                for ss=1:(2^level)*N1+1
                    idx=ss:NN1+1:NN1*NN2;
                    x_f=I_d(idx);
                    count(idx)=count(idx)+1;
                    [c,l]=wavedec(x_f,level,filter_fri);
                    y_n=c(1:l(1)).';
                    nn=length(y_n);
                    n_org=length(x_f);
                    kk=ceil(nn./4);

                    [t,lambda,c_m_n_mdf]=setup_filter(filter_fri,nn,level,a);
                    [t_org,~,~]=setup_filter(filter_fri,nn,level,a);
                    T_s=t(2)-t(1);
                    inloc_est=find_polyCons(kk,y_n,lambda,c_m_n_mdf)-loc_shift*T_s;
                    if up_flag
                        r_upres=round((inloc_est-t_org(1))/(t_org(2)-t_org(1))*upres+1);
                        r_upres=sort(r_upres((r_upres>1)&(r_upres<n_org*upres))); 
                        rr_est=unique(r_upres);
                        [Corg,Sorg]=wavedec(zeros(n_org,1),level,filter_fri);
                        [C_up,S_up]=wavedec(zeros(n_org*upres,1),log2(upres),filter_fri);
                        lin_rec_up=waverec([x_f(:).' zeros(1,length(C_up)-Sorg(1))],S_up,filter_fri);
                        x_est_up=linfit(rr_est,lin_rec_up);
                        x_rec_up=wavedec(x_est_up,log2(upres),filter_fri);
                        x_rec=x_rec_up(1:n_org);
                    else
                        r_est=round((inloc_est-t_org(1))/(t_org(2)-t_org(1))+1);
                        r_est=sort(r_est((r_est>1)&(r_est<NN)));
                        rr_est=unique(r_est);
                        x_rec=linfit(rr_est,x_f);
                    end 
                    Coef_d1(idx)=Coef_d1(idx)+x_rec';
                end
                Coef_d1=Coef_d1./count;

                %%%%%%% first downsampling and then upsampling line by line using FRI : diagonal line -45 degree %%%%%%%%
                Coef_d2=zeros(NN1,NN2);count=zeros(NN1,NN2);
                I_d=fliplr(I_f);
                for ss=1:(2^level)*N1+1
                    idx=ss:NN1+1:NN1*NN2;
                    x_f=I_d(idx);
                    count(idx)=count(idx)+1;
                    X_rec=[];
                    [c,l]=wavedec(x_f,level,filter_fri);
                    y_n=c(1:l(1)).';
                    nn=length(y_n);
                    n_org=length(x_f);
                    yn_org=y_n;
                    kk=ceil(nn./4);

                    [t,lambda,c_m_n_mdf]=setup_filter(filter_fri,nn,level,a);
                    [t_org,~,~]=setup_filter(filter_fri,nn,level,a);
                    T_s=t(2)-t(1);
                    inloc_est=find_polyCons(kk,y_n,lambda,c_m_n_mdf)-loc_shift*T_s;
                    if up_flag
                        r_upres=round((inloc_est-t_org(1))/(t_org(2)-t_org(1))*upres+1);
                        r_upres=sort(r_upres((r_upres>1)&(r_upres<n_org*upres)));
                        rr_est=unique(r_upres);
                        [Corg,Sorg]=wavedec(zeros(n_org,1),level,filter_fri);
                        [C_up,S_up]=wavedec(zeros(n_org*upres,1),log2(upres),filter_fri);
                        lin_rec_up=waverec([x_f(:).' zeros(1,length(C_up)-Sorg(1))],S_up,filter_fri);
                        x_est_up=linfit(rr_est,lin_rec_up);
                        x_rec_up=wavedec(x_est_up,log2(upres),filter_fri);
                        x_rec=x_rec_up(1:n_org);
                    else
                        r_est=round((inloc_est-t_org(1))/(t_org(2)-t_org(1))+1);
                        r_est=sort(r_est((r_est>1)&(r_est<NN)));
                        rr_est=unique(r_upres);
                        x_rec=linfit(rr_est,x_f);
                    end
                    Coef_d2(idx)=Coef_d2(idx)+x_rec';
                end
                Coef_d2=Coef_d2./count;
                Coef_d2=fliplr(Coef_d2);

                Coef_d1=range0toN(Coef_d1,rangeImg);
                Coef_d2=range0toN(Coef_d2,rangeImg);
                %%% select patches from either Coef_d1 or Coef_d2 for reconstruction according to the patch's gradient direction %%%
                [Gmag,Gdir] = imgradient(I_f);
                Coef_all=zeros(NN1,NN2);shift=1;psize=4;
                count=zeros(NN1,NN2);
                ang=[-135 45 -45 135];eggs=6;
                for xx=0:(NN1-psize)/shift
                    for yy=0:(NN2-psize)/shift
                        mag=Gmag(xx*shift+(1:psize),yy*shift+(1:psize));
                        dir=Gdir(xx*shift+(1:psize),yy*shift+(1:psize));
                        [value,idx]=sort(mag(:),'descend');
                        p_dir=(dir(idx(1:eggs)));id_vec=zeros(size(ang));
                        for mm=1:eggs
                            [~,id]=min(abs(p_dir(mm)-ang));
                            id_vec(id(1))=id_vec(id(1))+1;
                        end
                        [~,idm]=max(id_vec);
                        patch=zeros(NN1,NN2);

                        if (idm==1 || idm==2)
                            patch(xx*shift+(1:psize),yy*shift+(1:psize))=(Coef_d1(xx*shift+(1:psize),yy*shift+(1:psize)));
                        elseif (idm==3 || idm==4)
                            patch(xx*shift+(1:psize),yy*shift+(1:psize))=(Coef_d2(xx*shift+(1:psize),yy*shift+(1:psize)));
                        end
                        Coef_all=Coef_all+patch;
                        count(xx*shift+(1:psize),yy*shift+(1:psize))=count(xx*shift+(1:psize),yy*shift+(1:psize))+1;
                    end
                end
                Coef_all=Coef_all./count;
                Coef_all=range0toN(Coef_all,rangeImg);
                if (it==2)&&(supres==1)
                    level=2;
                else
                    level=1;
                end
                [c1,l1] = wavedec2(Coef_all,level,filter);
                [h1,v1,d1] = detcoef2('all',c1,l1,1);
                if (it==2)&&(supres==1)
                    [h2,v2,d2] = detcoef2('all',c1,l1,2);
                    c_rec=[4*I_int(:); h2(:); v2(:); d2(:); h1(:); v1(:); d1(:);];
                else
                    c_rec=[I_low(:);  h1(:); v1(:); d1(:);];
                end

                [c_org,l_org] = wavedec2(zeros(NN1,NN2),level,filter);
                I_rec = waverec2(c_rec,l_org,filter);
                I_rec=range0toN(I_rec,rangeImg);

                I_f=I_rec;
                I_diag=I_f;
            end

            if supres==0 
               I_fri{it}=I_rec;
            elseif supres==1
               I_fri2{it}=I_rec;
            else
               I_fri1{it}=I_rec;
            end

            if supres==1
                %%%%% CREATE INTERMEDIATE LEVEL %%%%%
                I_org0{it}=intp_bic(I_org{it},b);
                I_fri0{it}=I_fri{it};%intp_bic(I_fri{it},b);

                %%%%% CORRECT CURRENT FRI VERSION I_FRI2{IT} USING THE INTEMEDIATE LEVEL
                [I_up1{it},IDX,MSE]=FRI_errRepair_ann_lm(I_org0{it},I_fri0{it},I_fri2{it},it);
                I_up1{it}=range0toN(I_up1{it},rangeImg);

                level=it;

                [c1,l1] = wavedec2(I_up1{it},level,filter);
                [c_org,l_org] = wavedec2(zeros(NN1,NN2),level,filter);
                [h1,v1,d1] = detcoef2('all',c1,l1,1);

                if (it==2)
                [h2,v2,d2] = detcoef2('all',c1,l1,2);
                c_rec=[4*I_int(:); h2(:); v2(:); d2(:); h1(:); v1(:); d1(:);];
                else
                c_rec=[I_low(:);  h1(:); v1(:); d1(:);];

                end
                I_rec = waverec2(c_rec,l_org,filter);
                I_rec=range0toN(I_rec,rangeImg);
                I_f=I_rec; 
                %b=b+0.2;
                % I_tmp=intp_bic(I_rec,b);
                % [c1,l1] = wavedec2(I_tmp,2,filter);
                % I_up1to2{it}=appcoef2(c1,l1,filter,1)./2;
                I_up1to2{it}=intp_bic(I_rec,b/2);
            end
            if supres>=2
                %%%%%%%% UPDATE INTERMEDIATE LEVEL AND RE-CORRECT THE CUREENT FRI VERSION
                I_org1{it}=I_up1to2{it};
                %if (supres==2 || supres==4)
                I_up2{it}=FRI_errRepair_ann_lm2(I_org1{it},I_fri1{it},I_fri2{it},it,IDX,MSE);
                %else
                %[I_up2{it},IDX,MSE]=FRI_errRepair_ann_lm(I_org1{it},I_fri1{it},I_fri2{it},it);
                %end
                I_up2{it}=range0toN(I_up2{it},rangeImg);

                level=it;

                [c1,l1] = wavedec2(I_up2{it},level,filter);
                [c_org,l_org] = wavedec2(zeros(size(I_up2{it})),level,filter);
                [h1,v1,d1] = detcoef2('all',c1,l1,1);

                if (it==2)
                [h2,v2,d2] = detcoef2('all',c1,l1,2);

                c_rec=[4*I_int(:); h2(:); v2(:); d2(:); h1(:); v1(:); d1(:);];
                else
                c_rec=[2*I_org{it}(:);  h1(:); v1(:); d1(:);];

                end
                I_rec = waverec2(c_rec,l_org,filter);
                I_rec=range0toN(I_rec,rangeImg);
                I_f=I_rec;


                I_up1to2{it}=intp_bic(I_rec,b/2);
            end
        end
        I_upres{it}=I_rec;
    end
    if atf
        bd=4;
        PSNR_rec = PSNR(I_rec(bd+1:end-bd,bd+1:end-bd),I(bd+1:end-bd,bd+1:end-bd),255);
        disp(PSNR_rec);
        formatSpec = 'PSNR of %12s with filter %10s is %.2f dB\n';
        fileID = fopen('PSNR.txt','w');
        fprintf(fileID,formatSpec,imagename,filter_fri,PSNR_rec);
        fclose(fileID);
    end
    toc
    if layer == 3
        I_c = zeros(size(I_upres{it},1),size(I_upres{it},2),3);
        I_c(:,:,1)=uint8(I_upres{it});
        I_c(:,:,2)=uint8(hImcb);
        I_c(:,:,3)=uint8(hImcr);
        I_c=ycbcr2rgb(uint8(I_c));
        figure;
        imshow(I_c,[]);
        imwrite(uint8(I_c),'temp.bmp');
%         imwrite(uint8(I_c),'temp.jpg');
    else
        figure;
        imshow(I_rec,[]);
        imwrite(uint8(I_rec),[imagename '_PSNR' num2str(PSNR_rec,'%.2f') '_lm' num2str(lm_iter) '_diag' num2str(diag_iter) '_intbic_' filter '.bmp']);
%         imwrite(uint8(I_rec),[imagename '_PSNR' num2str(PSNR_rec,'%.2f') '_lm' num2str(lm_iter) '_diag' num2str(diag_iter) '_intbic_' filter '.jpg']);
    end
end

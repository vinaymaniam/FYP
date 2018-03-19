
%%
videos = [dir('data/*.mat');dir('data/*.y4m')];
numframes = 16;
downsample = 0;
crop_space = 0;
factorSR = 2;
factorI = 1.625;
lm_with_motion = 1;

do_bic_intp = 1;
do_temporal_sr = 1;
do_lm = 1;
do_updated_lm = 0;
do_color = 1;

do_compile = 1;

folder = [strrep(int2str(fix(clock)),' ',''),'_color',int2str(do_color),...
    '_crop',int2str(crop_space),'_downsmpl',int2str(downsample),...
    '_factor',int2str(factorSR),'_frames',int2str(numframes),...
    '_lm',int2str(do_lm),int2str(do_updated_lm),...
    '_lm_with_mot',int2str(lm_with_motion),'/'];
mkdir(folder);

%%
if do_compile
    if isunix
        coptimflags = '-O3 -ffast-math';
    else
        coptimflags = '/O2gtx /fp:fast';
    end
    eval(['mex COPTIMFLAGS="$COPTIMFLAGS ',coptimflags,'" ',...
        ' -DMEX_COMPILE_FLAG=1 -DWITH_JPEG=0 ',...
        '-output video_upsmpl_reg src/main_video_reg.cpp -lmwlapack -lmwblas' ]);
    eval(['mex COPTIMFLAGS="$COPTIMFLAGS ',coptimflags,'" ',...
        ' -DMEX_COMPILE_FLAG=1 -DWITH_JPEG=0 ',...
        '-output video_upsmpl_intp src/main_video_intp.cpp -lmwlapack -lmwblas' ]);
    eval(['mex COPTIMFLAGS="$COPTIMFLAGS ',coptimflags,'" ',...
        ' -DMEX_COMPILE_FLAG=1 -DWITH_JPEG=0 ',...
        '-output video_upsmpl_lm src/main_video_lm.cpp -lmwlapack -lmwblas' ]);
    clear coptimflags
end

%%
for video_idx=1:length(videos)
    %%
    disp(['Reading video ',videos(video_idx).name,'...']);
    filename_out = [folder,videos(video_idx).name(1:end-4),'_up.mat'];
    ext = videos(video_idx).name(end-2:end);
    if strcmpi(ext,'y4m')
        video = yuv4mpeg2mov(['data/',videos(video_idx).name]);
        disp('Converting to YUV space...');
        y = zeros(size(video(1).cdata,1),size(video(1).cdata,2),3,length(video));
        for i=1:length(video)
            y(:,:,:,i) = rgb2ycbcr(im2double(video(i).cdata));
        end
        crop_time = numframes*2;
        decimation = 1;
    elseif strcmpi(ext,'mat')
        load(['data/',videos(video_idx).name],'y')
        y = im2double(y);
        for i=1:size(y,4)
            y(:,:,:,i) = rgb2ycbcr(y(:,:,:,i));
        end
        crop_time = numframes;
        decimation = 0;
    else
        error('extension not supported')
    end
    
    %%
    if downsample
        disp('Downsampling video');
        dwtmode('per','nodisp');
        tmp = [];
        for i=1:size(y,4)
            for c=1:3
                frame = y(:,:,c,i);
                for n=1:downsample
                    frame = dwt2(frame,'haar')/2;
                end
                tmp(:,:,c,i) = frame;
            end
        end
        y = tmp;
    end

    %%
    if crop_time
        disp('Cropping video along time...');
        t0 = 1;
        y = y(:,:,:,t0:min(size(y,4),t0+crop_time-1));
    end
    if crop_space
        disp('Cropping video in space...');
        if crop_space==1
            crop_space = 2;
        end
        mr = round(size(y,1)/crop_space);
        mc = round(size(y,2)/crop_space);
        mmr = round((size(y,1)-mr)/2);
        mmc = round((size(y,2)-mc)/2);
        y = y(mmr:mmr+mr-1,mmc:mmc+mc-1,:,:);
    end
    if decimation
        disp('Decimating video...');
        z = y(:,:,:,1:factorSR:end);
    else
        z = y;
    end
    save(filename_out, 'z', 'y')
    
    %%
    if do_color
        y_rgb = zeros(size(y));
        for i=1:size(y,4)
            y_rgb(:,:,:,i) = ycbcr2rgb(y(:,:,:,i));
        end
        save(filename_out, 'y_rgb', '-append')
    end
    
    %%
    if do_bic_intp
        disp('Bicubic interpolation video...');
        y_bic_yuv = zeros(size(z).*[1 1 1 factorSR]);
        x = 0:size(z,4)-1;
        xv = (0:factorSR*size(z,4)-1)/factorSR;
        for i=1:size(z,1)
            for j=1:size(z,2)
                for c=1:3
                    v = squeeze(z(i,j,c,:));
                    xsec = interp1(x,v, xv, 'pchip');
                    xsec(xv>x(end)) = v(end);
                    y_bic_yuv(i,j,c,:) = xsec;
                end
            end
        end
        y_bic = squeeze(y_bic_yuv(:,:,1,:));
        save(filename_out, 'y_bic', '-append')
        if do_color
            y_bic_rgb = zeros(size(y_bic_yuv));
            for i=1:size(y_bic_rgb,4)
                y_bic_rgb(:,:,:,i) = ycbcr2rgb(cat(3, y_bic_yuv(:,:,1,i),y_bic_yuv(:,:,2,i),y_bic_yuv(:,:,3,i)));
            end
            save(filename_out, 'y_bic_rgb', '-append')
        end
    end

    %%

    if do_temporal_sr
        if do_color
            disp(['Temporal upsampling color video [',num2str(size(z)),']...']);
            tic;
            disp('Luminance Y')
            [MF,y_reg,y_intp] = video_upsmpl_reg(squeeze(z(:,:,1,:)), factorSR);
            disp('Chrominance Cb')
            [y_intp_Cb] = video_upsmpl_intp(MF, squeeze(z(:,:,2,:)), factorSR);
            disp('Chrominance Cr')
            [y_intp_Cr] = video_upsmpl_intp(MF, squeeze(z(:,:,3,:)), factorSR);
            toc
            y_reg_rgb = zeros(size(y_reg,1),size(y_reg,2),3,size(y_reg,3));
            y_intp_rgb = zeros(size(y_reg_rgb));
            for i=1:size(y_reg_rgb,4)
                y_reg_rgb(:,:,:,i) = ycbcr2rgb(cat(3, y_reg(:,:,i),y_intp_Cb(:,:,i),y_intp_Cr(:,:,i)));
                y_intp_rgb(:,:,:,i) = ycbcr2rgb(cat(3, y_intp(:,:,i),y_intp_Cb(:,:,i),y_intp_Cr(:,:,i)));
            end
            save(filename_out, 'MF', 'y_reg', 'y_intp', ...
                'y_intp_Cb', 'y_intp_Cr', 'y_reg_rgb', 'y_intp_rgb', '-append')
        else
            disp(['Temporal upsampling grayscale video [',num2str(size(z)),']...']);
            tic;
            [MF,y_reg,y_intp] = video_upsmpl_reg(squeeze(z(:,:,1,:)), factorSR);
            toc
            save(filename_out, 'MF', 'y_reg', 'y_intp', '-append')
        end
    end

    %%
    if do_temporal_sr && do_lm
        disp(['Linear Mapping video [',num2str(size(y_reg)),']...']);
        tic;
        [y_lm,y_inter,y_inter_est,y_inter_intp] = video_upsmpl_lm(...
            squeeze(z(:,:,1,:)),y_reg,factorI,lm_with_motion);
        toc
        save(filename_out, 'y_lm', 'y_inter','y_inter_est', '-append')
        if do_color
            y_lm_rgb = zeros(size(y_lm,1),size(y_lm,2),3,size(y_lm,3));
            for i=1:size(y_lm_rgb,4)
                y_lm_rgb(:,:,:,i) = ycbcr2rgb(cat(3, y_lm(:,:,i),y_intp_Cb(:,:,i),y_intp_Cr(:,:,i)));
            end
            save(filename_out, 'y_lm_rgb', '-append')
        end
    end
    
    %%
    if do_temporal_sr && do_lm && do_updated_lm
        disp(['Updated Linear Mapping video [',num2str(size(y_reg)),']...']);
        tic;
        [y_lm_hat, y_inter_hat, y_inter_est_hat, y_inter_intp_hat] = ...
            video_upsmpl_lm(y_lm,y_reg,factorI,lm_with_motion);
        toc
        save(filename_out, 'y_lm_hat','y_inter_hat', 'y_inter_est_hat', '-append')
        if do_color
            y_lm_hat_rgb = zeros(size(y_lm_hat,1),size(y_lm_hat,2),3,size(y_lm_hat,3));
            for i=1:size(y_lm_hat_rgb,4)
                y_lm_hat_rgb(:,:,:,i) = ycbcr2rgb(cat(3, y_lm_hat(:,:,i),y_intp_Cb(:,:,i),y_intp_Cr(:,:,i)));
            end
            save(filename_out, 'y_lm_hat_rgb', '-append')
        end
    end
    disp(' ')

end


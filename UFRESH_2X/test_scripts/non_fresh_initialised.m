clear;
dwtmode('per')

addpath('vinay')
addpath('../Python/data_files')

%% Metric being tested here
% nvals = [128,256,512,1024,2048,4096,8192,16384];
nvals = 16384;
% ---------------------------------------------------------------------------
meanpsnrs = zeros(length(nvals),2);
meanssims = zeros(length(nvals),2);
meantimeperpixel = zeros(length(nvals),2);
for set = 1:2
    pattern = '*.bmp';
    directory_y = sprintf('Testing_Images/GT/Set%i',9*set-4);     
    YpathCell = glob(directory_y, pattern );
    Ycell = load_images( YpathCell );
    for i = 1:length(Ycell)
        tbic = tic;
        y = Ycell{i};
        x1 = imresize(imresize(y,0.5),2);
        x = iterativeBackproj(x1);
        mask = ones(size(y)) - [zeros(3,size(y,2));...
            [zeros(size(y,1)-5,3),ones(size(y,1)-5,size(y,2)-5),zeros(size(y,1)-5,2)];...
            zeros(2,size(y,2))];
        mask = logical(mask);
        x(mask) = x1(mask);
        psnr(x,y)
        Xcell{i} = x;
        toc(tbic)/numel(Xcell{i});
    end


    blocksize = [5, 5]; % the size of patch.
    stepsize = [1, 1];  
    if length(Xcell) ~= length(Ycell)	
        error('Error: The number of X images is not equal to the number of Y images!');
    end
    for n = nvals
        fprintf('################   %d    #####################\n',n)
        postpsnr=zeros(1,length(Xcell)); prepsnr = zeros(1,length(Xcell));
        postssim=zeros(1,length(Xcell)); pressim = zeros(1,length(Xcell));
        tpp = zeros(1,length(Xcell));
        %% Begin SR
        for imgIdx = 1:length(Xcell)
            stopwatch1 = tic;
            fprintf('--------------------------------------------------------\n')
            fprintf('Processing image %d of total %d ... \n', imgIdx, length(Xcell));

            Xtest = Xcell{imgIdx}; % LowRresolution image X            
            Ytest = Ycell{imgIdx}; % HighResolution image Y
            Xtest = imresize(Xtest,size(Ytest));
            fprintf('[BEFORE] PSNR = %.1f     SSIM = %.3f\n', psnr(Xtest,Ytest),ssim(Xtest,Ytest));
            prepsnr(imgIdx) = psnr(Xtest,Ytest);
            pressim(imgIdx) = ssim(Xtest,Ytest);    
            %% NEXT STEP - IMPLEMENT RESIDUAL LEARNING(FROM FRESH) IN V2
            %% NEED TO FIGURE OUT A WAY TO USE Res IN pyMapCell TO DO SOME SORT OF ERROR CORRECTION
            for stage = 1:1
                %% Load trained model
                load(sprintf('%ipyHeirarchy%i_NF',stage,n));
                heirarchy = single(heirarchy);   
                load(sprintf('%ipyMap%icell192_NF',stage,n));
                ensembleSize = 4; % low ensemble size --> not too big of a drop in quality
                Xrec = zeros([size(Xtest),ensembleSize]);
                for rot = 1:ensembleSize
                    X = affine(Xtest,abs(mod(rot-1,4)),rot-4,0);
                    X = ufresh2(X, blocksize, heirarchy, index, Map);
                    X = affine(X,abs(mod(rot-1,4)),rot-4,1);
                    X = range0toN(X,[0,1]);
                    Xrec(:,:,rot) = X;            
                end        
                Xtest = mean(Xrec,3);
                % Clip image to 0-1 range
                Xtest = range0toN(Xtest,[0,1]);
            end
            fprintf('[AFTER]  PSNR = %.1f     SSIM = %.3f\n', psnr(Xtest,Ytest),ssim(Xtest,Ytest));         
            postpsnr(imgIdx)=psnr(Xtest,Ytest); 
            postssim(imgIdx)=ssim(Xtest,Ytest); 
            toc(stopwatch1)
            tpp(imgIdx) = toc(stopwatch1)/numel(Ytest);
    %         visualiseImprovement(Xcell{imgIdx}, Xtest, Ytest)
        end
        fprintf('============================================================\n')
        fprintf('Average PSNR across all images = %.2f\n',mean(postpsnr))
        fprintf('Average SSIM across all images = %.3f\n',mean(postssim))
        fprintf('Average improvement in PSNR = %.2f\n',mean(postpsnr-prepsnr))
        fprintf('Average improvement in SSIM = %.3f\n',mean(postssim-pressim))
        fprintf('Average Time Per Pixel = %f\n',mean(tpp))
        meanpsnrs(nvals == n,set) = mean(postpsnr);
        meanssims(nvals == n,set) = mean(postssim);
        meantimeperpixel(nvals == n,set) = mean(tpp);    
    end

    plt = 0;
    if plt
        meanpsnrs = (5*meanpsnrs(:,1)+14*meanpsnrs(:,2))/19;
        meanssims = (5*meanssims(:,1)+14*meanssims(:,2))/19;
        meantimeperpixel = (5*meantimeperpixel(:,1)+14*meantimeperpixel(:,2))/19;
        % Plot results
        figure;
        plot(nvals, meanpsnrs)
        xlabel('Number of Centroids')
        ylabel('Average PSNR(dB)')
        title('Impact of # Centroids on PSNR(Set 5)');
        % title('Impact of # Centroids on PSNR(Set 14)');

        figure;
        plot(nvals, meanssims)
        xlabel('Number of Centroids')
        ylabel('Average SSIM')
        title('Impact of # Centroids on SSIM(Set 5)');
        % title('Impact of # Centroids on SSIM(Set 14)');

        figure;
        plot(nvals, meantimeperpixel,'o--')
        xlabel('Number of Centroids')
        ylabel('Average Time Per Pixel')
        title('Impact of # Centroids on Runtime Speed');
    end
end


        
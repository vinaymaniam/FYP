clear;
dwtmode('per')

addpath('vinay')
addpath('../Python/data_files')

%% Metric being tested here
nvals = [1,2,3,4];
% ---------------------------------------------------------------------------
meanpsnrs = zeros(length(nvals),2);
meanssims = zeros(length(nvals),2);
meantimeperpixel = zeros(length(nvals),2);

n = 8192;
for set = 1:2
    directory_x = sprintf('Testing_Images/FRESH_upscaled/Set%i',9*set-4);     
    directory_y = sprintf('Testing_Images/GT/Set%i',9*set-4); 
    
    pattern = '*.bmp';
    XpathCell = glob(directory_x, pattern );
    Xcell = load_images( XpathCell );
    YpathCell = glob(directory_y, pattern );
    Ycell = load_images( YpathCell );
    fprintf('=================== SET %i ===========================\n',9*set-4)
    blocksize = [5, 5]; % the size of patch.
    stepsize = [1, 1];  
    if length(Xcell) ~= length(Ycell)	
        error('Error: The number of X images is not equal to the number of Y images!');
    end
    
    postpsnr=zeros(4,length(Xcell)); prepsnr = zeros(4,length(Xcell));
    postssim=zeros(4,length(Xcell)); pressim = zeros(4,length(Xcell));
    tpp = zeros(1,length(Xcell));
    %% Begin SR
    for imgIdx = 1:length(Xcell)
        stopwatch1 = tic;
        fprintf('--------------------------------------------------------\n')
        fprintf('Processing image %d of total %d ... \n', imgIdx, length(Xcell));

        Xtest = Xcell{imgIdx}; % LowRresolution image X
        Ytest = Ycell{imgIdx}; % HighResolution image Y    
        fprintf('[BEFORE] PSNR = %.1f     SSIM = %.3f\n', psnr(Xtest,Ytest),ssim(Xtest,Ytest));             
        %% NEXT STEP - IMPLEMENT RESIDUAL LEARNING(FROM FRESH) IN V2
        %% NEED TO FIGURE OUT A WAY TO USE Res IN pyMapCell TO DO SOME SORT OF ERROR CORRECTION
        for stage = 1:3
            prepsnr(stage,imgIdx) = psnr(Xtest,Ytest);
            pressim(stage,imgIdx) = ssim(Xtest,Ytest);   
            %% Load trained model
            load(sprintf('%ipyHeirarchy%i',stage,n));
            heirarchy = single(heirarchy);   
            load(sprintf('%ipyMap%icell192',stage,n));
            ensembleSize = 4; % low ensemble size --> not too big of a drop in quality
            Xrec = zeros([size(Xtest),ensembleSize]);
            for rot = 1:ensembleSize
                X = rot90(Xtest, (rot-1));                        
                X = ufresh2(X, blocksize, heirarchy, index, Map);
                X = rot90(X, 4-(rot-1));
                X = range0toN(X,[0,1]);
                Xrec(:,:,rot) = X;            
            end        
            Xtest = mean(Xrec,3);
            Xtest = backprojection_2X(Xtest, Ytest, 'bior4.4'); 
            % Clip image to 0-1 range
            Xtest = range0toN(Xtest,[0,1]);
            %%
            fprintf('[AFTER]  PSNR = %.1f     SSIM = %.3f\n', psnr(Xtest,Ytest),ssim(Xtest,Ytest));         
            postpsnr(stage,imgIdx)=psnr(Xtest,Ytest); 
            postssim(stage,imgIdx)=ssim(Xtest,Ytest); 
        end            
        toc(stopwatch1)
        tpp(imgIdx) = toc(stopwatch1)/numel(Ytest);
%         visualiseImprovement(Xcell{imgIdx}, Xtest, Ytest)
    end
    fprintf('============================================================\n')        
    meanpsnrs(:,set) = mean(postpsnr,2);
    meanssims(:,set) = mean(postssim,2);
    meantimeperpixel(:,set) = mean(tpp);    

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


        
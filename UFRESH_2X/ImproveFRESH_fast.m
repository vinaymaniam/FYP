clear;
dwtmode('spd')

addpath('vinay')
% testimgs = 'bior44';
testimgs = 'fresh';
        
directory_x = 'Testing_Images/FRESH_upscaled/Set5'; 
% directory_x = 'Testing_Images/bicubic/Set5'; 
% directory_x = 'Testing_Images/bior44/Set5';
% directory_x = 'Testing_Images/db2/Set5';
pattern = '*.bmp';
directory_y = 'Testing_Images/GT/Set5'; 

XpathCell = glob(directory_x, pattern );
Xcell = load_images( XpathCell );
YpathCell = glob(directory_y, pattern );
Ycell = load_images( YpathCell );

switch testimgs
    case 'bior44'
        zcell = cell(length(Xcell));
        for i = 1:length(Xcell)
            zcell{i} = Xcell{i};
            Xcell{i} = imresize(Xcell{i}, size(Ycell{i}));
        end
end

blocksize = [5, 5]; % the size of patch.
% If stepsize < blocksize, there exists overlapping aeras among adjacent patches.
stepsize = [1, 1];  
if length(Xcell) ~= length(Ycell)	
	error('Error: The number of X images is not equal to the number of Y images!');
end
postpsnr=zeros(1,length(Xcell)); prepsnr = zeros(1,length(Xcell));
postssim=zeros(1,length(Xcell)); pressim = zeros(1,length(Xcell));
for imgIdx = 1:length(Xcell)
    stopwatch1 = tic;
	fprintf('--------------------------------------------------------\n')
	fprintf('Processing image %d of total %d ... \n', imgIdx, length(Xcell));

    Xtest = Xcell{imgIdx}; % LowRresolution image X
    Ytest = Ycell{imgIdx}; % HighResolution image Y    
    fprintf('[BEFORE] PSNR = %.1f     SSIM = %.3f\n', psnr(Xtest,Ytest),ssim(Xtest,Ytest));
    prepsnr(imgIdx) = psnr(Xtest,Ytest);
    pressim(imgIdx) = ssim(Xtest,Ytest);
    %% Load trained model
    load(sprintf('pyHeirarchy4096'));
    heirarchy = single(heirarchy);   
    % python uses 0 indexing
    % index = index + 1;        
    load(sprintf('pyMap4096cell96'));    
    filt = 'db2'; % db2 gives much better results for FRESH input
    % Cascading makes things WORSE
    %% NEXT STEP - IMPLEMENT RESIDUAL LEARNING(FROM FRESH) IN V2
    for stage = 1:1%2
        ensembleSize = 4; % low ensemble size --> not too big of a drop in quality
        Xrec = zeros([size(Xtest),ensembleSize]);
        for rot = 1:ensembleSize
            X = rot90(Xtest, (rot-1));                        
            X = ufresh2(X, blocksize, heirarchy, index, Map);
            X = rot90(X, 4-(rot-1));
            X = backprojection_2X(X, Ytest, filt);
            X = range0toN(X,[0,1]);
            Xrec(:,:,rot) = X;            
        end        
%         res = resLearn(Xtest, 2, blocksize, heirarchy, index, Map, filt);
        Xtest = mean(Xrec,3);
%         Xtest = Xtest + res;
        Xtest = backprojection_2X(Xtest, Ytest, filt);                
        Xtest = range0toN(Xtest,[0,1]);
    end
    fprintf('[AFTER]  PSNR = %.1f     SSIM = %.3f\n', psnr(Xtest,Ytest),ssim(Xtest,Ytest));
    postpsnr(imgIdx)=psnr(Xtest,Ytest); 
    postssim(imgIdx)=ssim(Xtest,Ytest); 
    toc(stopwatch1)
end
fprintf('============================================================\n')
fprintf('Average PSNR across all images = %.2f\n',mean(postpsnr))
fprintf('Average improvement in PSNR = %.2f\n',mean(postpsnr-prepsnr))
fprintf('Average improvement in SSIM = %.2f\n',mean(postssim-pressim))
        
clear;
dwtmode('spd')

addpath('vinay')
addpath('../Python/data_files')
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
stepsize = [1, 1];  
if length(Xcell) ~= length(Ycell)	
	error('Error: The number of X images is not equal to the number of Y images!');
end
postpsnr=zeros(1,length(Xcell)); prepsnr = zeros(1,length(Xcell));
postssim=zeros(1,length(Xcell)); pressim = zeros(1,length(Xcell));
%% Load trained model
n = 1024;
load(sprintf('pyHeirarchy%i',n));
heirarchy = single(heirarchy);   
load(sprintf('pyMap%icell96',n));    
filt = 'db2'; % db2 gives much better results for FRESH input
%% Begin SR
for imgIdx = 1:length(Xcell)
    stopwatch1 = tic;
	fprintf('--------------------------------------------------------\n')
	fprintf('Processing image %d of total %d ... \n', imgIdx, length(Xcell));

    Xtest = Xcell{imgIdx}; % LowRresolution image X
    Ytest = Ycell{imgIdx}; % HighResolution image Y    
    fprintf('[BEFORE] PSNR = %.1f     SSIM = %.3f\n', psnr(Xtest,Ytest),ssim(Xtest,Ytest));
    prepsnr(imgIdx) = psnr(Xtest,Ytest);
    pressim(imgIdx) = ssim(Xtest,Ytest);    
    %% NEXT STEP - IMPLEMENT RESIDUAL LEARNING(FROM FRESH) IN V2
    %% NEED TO FIGURE OUT A WAY TO USE Res IN pyMapCell TO DO SOME SORT OF ERROR CORRECTION
    for stage = 1:1%2 %Cascading stages actually makes things worse unless you train separate model for each stage
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
        Xtest = backprojection_2X(Xtest, Ytest, filt);     
        % Clip image to 0-1 range
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
        
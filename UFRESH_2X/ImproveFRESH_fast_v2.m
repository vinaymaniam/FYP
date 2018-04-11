clear;
dwtmode('sym')
% addpath('ksvd');
% addpath('ksvd/ksvdbox');
% addpath('ksvd/ksvdbox/private_ccode');
% addpath('ksvd/ompbox');
addpath('vinay')

% directory_x = 'Testing_Images/FRESH_upscaled/Set5'; 
% directory_x = 'Testing_Images/bicubic/Set5';
directory_x = 'Testing_Images/bior44/Set5';
% directory_x = 'Testing_Images/db2/Set5';
pattern = '*.bmp';
directory_y = 'Testing_Images/GT/Set5'; 

XpathCell = glob(directory_x, pattern );
Xcell = load_images( XpathCell );
YpathCell = glob(directory_y, pattern );
Ycell = load_images( YpathCell );

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
	fprintf('--------------     %d of %d     --------------\n', imgIdx, length(Xcell))
    
    Ytest = Ycell{imgIdx}; % HighResolution image Y 
    Xtest = imresize(Xcell{imgIdx}, size(Ytest)); % LowRresolution image X
    Ilowc = Xcell{imgIdx};
    
    fprintf('[BEFORE] PSNR = %.1f     SSIM = %.3f\n', psnr(Xtest,Ytest),ssim(Xtest,Ytest));
    prepsnr(imgIdx) = psnr(Xtest,Ytest);
    pressim(imgIdx) = ssim(Xtest,Ytest);
    %% Load trained model
    load(sprintf('pyHeirarchy4096'));
    heirarchy = single(heirarchy);
    load(sprintf('pyMap4096cell96'));   
    filt = 'bior4.4';
    for stage = 1:2
        ensembleSize = 4;
        Xrec = zeros([size(Xtest),ensembleSize]);
        for rot = 1:ensembleSize
            X = rot90(Xtest, rot-1);
            X = ufresh2(X,blocksize,heirarchy,index, Map);   
            % PROBLEM: ROTATING Ilowc NOT SAME AS GETTING Ilowc FROM
            % ROTATED Ytest
%             X = backprojection(X,rot90(Ilowc, rot-1), filt);
            X = rot90(X, 4-(rot-1));
            X = backprojection(X, Ilowc, filt);
            Xrec(:,:,rot) = X;
        end
        Xtest = mean(Xrec,3);
        Xtest = backprojection(Xtest, Ilowc, filt);        
    end
    fprintf('[AFTER]  PSNR = %.1f     SSIM = %.3f\n', psnr(Xtest,Ytest),ssim(Xtest,Ytest));
    postpsnr(imgIdx)=psnr(Xtest,Ytest); 
    postssim(imgIdx)=ssim(Xtest,Ytest); 
    toc(stopwatch1)
end
fprintf('============================================\n')
fprintf('Average PSNR across all images = %.1f\n',mean(postpsnr))
fprintf('Average SSIM across all images = %.3f\n',mean(postssim))
fprintf('Average improvement in PSNR = %.2f(%.1f%%)\n',mean(postpsnr-prepsnr), 100*mean(postpsnr-prepsnr)/mean(prepsnr))
fprintf('Average improvement in SSIM = %.2f(%.1f%%)\n',mean(postssim-pressim), 100*mean(postssim-pressim)/mean(pressim))
        
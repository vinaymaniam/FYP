clear;
% addpath('ksvd');
% addpath('ksvd/ksvdbox');
% addpath('ksvd/ksvdbox/private_ccode');
% addpath('ksvd/ompbox');
addpath('vinay')

directory_x = 'Testing_Images/FRESH_upscaled/Set5'; 
pattern = '*.bmp';
directory_y = 'Testing_Images/GT/Set5'; 

XpathCell = glob(directory_x, pattern );
Xcell = load_images( XpathCell );
YpathCell = glob(directory_y, pattern );
Ycell = load_images( YpathCell );

blocksize = [5, 5]; % the size of patch.
stepsize = [1, 1];  % the step of extracting image patches. If stepsize < stepsize, there exists overlapping aeras among adjacent patches.
if length(Xcell) ~= length(Ycell)	
	error('Error: The number of X images is not equal to the number of Y images!');
end
Psnr=zeros(1,length(Xcell)); prepsnr = zeros(1,length(Xcell));
for imgIdx = 1:length(Xcell)
    stopwatch1 = tic;
	fprintf('--------------------------------------------------------\n')
	fprintf('Processing image %d of total %d ... \n', imgIdx, length(Xcell));

    Xtest = Xcell{imgIdx}; % LowRresolution image X
    Ytest = Ycell{imgIdx}; % HighResolution image Y    
    fprintf('PSNR before processing = %.1f\n', psnr(Xtest,Ytest))
    prepsnr(imgIdx) = psnr(Xtest,Ytest);
    %% Load trained model
    load(sprintf('pyHeirarchy4096'));
    heirarchy = single(heirarchy);   
    % python uses 0 indexing
    % index = index + 1;        
    load(sprintf('pyMap4096cell96'));    
    for stage = 1:2
        Xrec = zeros([size(Xtest),4]);
        for rot = 1:4   
            Xtest_rot = imrotate(Xtest, 90*(rot-1));
            X = ufresh2(Xtest_rot,blocksize,heirarchy,index, Map);
            X = backprojection_2X(X,imrotate(Ytest, 90*(rot-1)),'db2');
            X = imrotate(X, 360-90*(rot-1));
            Xrec(:,:,rot) = X;            
        end
        Xtest = mean(Xrec,3);
        Xtest = backprojection_2X(Xtest, Ytest,'bior3.3');
    end
    fprintf('PSNR after processing = %.1f\n', psnr(Xtest,Ytest))
    Psnr(imgIdx)=psnr(Xtest,Ytest); 
    toc(stopwatch1)
end
fprintf('============================================================\n')
fprintf('Average PSNR across all images = %.1f\n',mean(Psnr))
fprintf('Average improvement in PSNR = %.1f\n',mean(Psnr-prepsnr))
        
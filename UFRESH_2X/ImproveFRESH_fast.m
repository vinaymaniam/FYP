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
Psnr=zeros(1,length(Xcell));
for imgIdx = 1:1%length(Xcell)
    stopwatch1 = tic;
	fprintf('--------------------------------------------------------\n')
	fprintf('Processing image %d of total %d ... \n', imgIdx, length(Xcell));

    Xtest = Xcell{imgIdx}; % LowRresolution image X
    Ytest = Ycell{imgIdx}; % HighResolution image Y    
    
    for stage = 1:2
        % Using custom centroids for each stage barely makes a difference
        % approximately 0.07dB on average(and this comes at the cost of 
        % DOUBLED training time to get the second codebook
        %% Now that training is so fast might be worth it to make custom
        % codebooks AND Maps for each stage as before
        load(sprintf('pyHeirarchy4096'));
        heirarchy = single(heirarchy);   
        % python uses 0 indexing
%         index = index + 1;
        
        load(sprintf('pyMap4096cell96'));
        Xrec = zeros([size(Xtest),4]);
        for rot = 1:4   
            Xtest_rot = imrotate(Xtest, 90*(rot-1));
            X = ufresh2(Xtest_rot,blocksize,stepsize,heirarchy,index, Map);
            X = imrotate(X, 360-90*(rot-1));
            X = backprojection_2X(X,Ytest);
            Xrec(:,:,rot) = X;            
        end
        ensembleMean = mean(Xrec,3);
        
        Xtest = backprojection_2X(ensembleMean, Ytest);
   
        clear Center Map
    end
%     psnr(Xtest(3:end-2,3:end-2),Ytest(3:end-2,3:end-2)) % PSNR calculation
    psnr(Xtest,Ytest) % PSNR calculation
%     Psnr(imgIdx)=psnr(Xtest(3:end-2,3:end-2),Ytest(3:end-2,3:end-2)); 
    Psnr(imgIdx)=psnr(Xtest,Ytest); 
    toc(stopwatch1)
end
mean(Psnr)
        
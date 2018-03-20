clear;
addpath('ksvd');
addpath('ksvd/ksvdbox');
addpath('ksvd/ksvdbox/private_ccode');
addpath('ksvd/ompbox');
addpath('utils')
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
for img_index = 1:length(Xcell)
    stopwatch1 = tic;
	disp('--------------------------------------------------------')
	fprintf('Processing image %d of total %d ... \n', img_index, length(Xcell));

    Xtest = Xcell{img_index}; % LowRresolution image X
    Ytest = Ycell{img_index}; % HighResolution image Y    
    
    for stage = 1:2
        % Using custom centroids for each stage barely makes a difference
        % approximately 0.07dB on average(and this comes at the cost of 
        % DOUBLED training time to get the second codebook
        % load(sprintf('Center2048_%istage',stage)); 
        % Center=single(Center');
        % load(sprintf('Center4096'));
        % Center = single(Center);
        load(sprintf('Heirarchy4096'));
        heirarchy = single(heirarchy);        
        
        %load(sprintf('Map2048_%istage',stage));
        load(sprintf('Map4096cell'));
        Xrec = zeros([size(Xtest),4]);
        for rot = 1:4   
            Xtest_rot = imrotate(Xtest, 90*(rot-1));
            % X = ufresh2(Xtest_rot,blocksize,stepsize,Center, Map);
            X = ufresh2(Xtest_rot,blocksize,stepsize,heirarchy,index, Map);
            X = imrotate(X, 360-90*(rot-1));
            X = backprojection_2X(X,Ytest);
            Xrec(:,:,rot) = X;            
        end
        ensembleMean = mean(Xrec,3);
        
        Xtest = backprojection_2X(ensembleMean, Ytest);
   
        clear Center Map
    end
    psnr(Xtest(3:end-2,3:end-2),Ytest(3:end-2,3:end-2)) % PSNR calculation
    Psnr(img_index)=psnr(Xtest(3:end-2,3:end-2),Ytest(3:end-2,3:end-2)); 
    toc(stopwatch1)
end
mean(Psnr)
        
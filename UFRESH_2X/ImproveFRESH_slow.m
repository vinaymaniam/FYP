clear;
addpath('ksvd');
addpath('ksvd/ksvdbox');
addpath('ksvd/ksvdbox/private_ccode');
addpath('ksvd/ompbox');
addpath('utils')

directory_x = 'Testing_Images/FRESH_upscaled/Set5'; 
pattern = '*.bmp';
directory_y = 'Testing_Images/GT/Set5'; 

XpathCell = glob(directory_x, pattern );
Xcell = load_images( XpathCell );

YpathCell = glob(directory_y, pattern );
Ycell = load_images( YpathCell );
ensemble=0; % 1: enable ensemble learning, 0: disable ensemble learning
blocksize = [5, 5]; % the size of patch.
stepsize = [1,1];  % the step of extracting image patches. If stepsize < stepsize, there exists overlapping aeras among adjacent patches.
if length(Xcell) ~= length(Ycell)	
	error('Error: The number of X images is not equal to the number of Y images!');
end
Psnr=zeros(1,length(Xcell));
for img_index = 1: length(Xcell)
    
	disp('--------------------------------------------------------')
	fprintf('Processing image %d of total %d ... \n', img_index, length(Xcell));
    
    load Center2048_1stage; Center=single(Center');
    load Map2048_1stage;
            
    X_test = Xcell{img_index}; % LR image X
	Y_test = Ycell{img_index}; % HR image Y
    if (ensemble) % whether using ensemble learning or not; enabled, can be very slow!
        X_test_90=imrotate(X_test,90);
        X_test_180=imrotate(X_test,180);
        X_test_270=imrotate(X_test,270);

        X_rec_im0=ufresh_slow(X_test,blocksize,stepsize,Center, Map);
        X_rec_im90=ufresh_slow(X_test_90,blocksize,stepsize,Center, Map);
        X_rec_im180=ufresh_slow(X_test_180,blocksize,stepsize,Center, Map);
        X_rec_im270=ufresh_slow(X_test_270,blocksize,stepsize,Center, Map);

        X_rec_im90=imrotate(X_rec_im90,270);
        X_rec_im180=imrotate(X_rec_im180,180);
        X_rec_im270=imrotate(X_rec_im270,90);

        X_rec_im0=backpropogation_2X(X_rec_im0,Y_test);%backprojection
        X_rec_im90=backpropogation_2X(X_rec_im90,Y_test);
        X_rec_im180=backpropogation_2X(X_rec_im180,Y_test);
        X_rec_im270=backpropogation_2X(X_rec_im270,Y_test);

        X_rec_im=( X_rec_im0+X_rec_im90+X_rec_im180+X_rec_im270)./4; % average the four HR images
        I_bp=backpropogation_2X(X_rec_im,Y_test); % the reconstructed HR image
    else
        X_rec_im=ufresh_slow(X_test,blocksize,stepsize,Center, Map);   
        I_bpbackpropogation_2X(X_rec_im,Y_test); % the reconstructed HR image
    end
     
     Psnr(I_bp(3:end-2,3:end-2),Y_test(3:end-2,3:end-2)) % PSNR calculation
     Psnr(img_index)=psnr(I_bp(3:end-2,3:end-2),Y_test(3:end-2,3:end-2));
end
mean(Psnr)
        
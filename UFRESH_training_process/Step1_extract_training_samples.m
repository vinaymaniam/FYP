clear;
addpath('ksvd');
addpath('ksvd/ksvdbox');
addpath('ksvd/ksvdbox/private_ccode');
addpath('ksvd/ompbox');
addpath('utils')

directory_x = 'TrainingData/FRESH_1stage'; 
pattern = '*.bmp';
directory_y = 'TrainingData/GT'; 

XpathCell = glob(directory_x, pattern );
Xcell = load_images( XpathCell );

YpathCell = glob(directory_y, pattern );
Ycell = load_images( YpathCell );
blocksize = [5, 5]; % the size of each image patch.
trainnum = 20000;
variance_Thresh = 0.02;
X = [];
Y = [];

[n1, n2] = size(Xcell);
	for i = 1: n1
		for j = 1: n2
			X_temp = Xcell{i,j};         
			Y_temp = Ycell{i,j};

			p = ndims(Y_temp); % Number of dimensions.
			ids = cell(p,1); % indices
			[ids{:}] = reggrid(size(Y_temp)-blocksize+1, trainnum, 'eqdist');
			Y_a = sampgrid(Y_temp,blocksize,ids{:});
			X_a = sampgrid(X_temp,blocksize,ids{:});
           
			X = [X X_a];            
			Y = [Y Y_a];
		end
	end
	% remove the mean (the dc component from each patch)
    Y = Y - repmat(mean(X), size(X, 1), 1);
	X = X - repmat(mean(X), size(X, 1), 1); % remove the low resolution 
	
    
    % Discard those patch pairs with too small variance.
	Xnorm2 = sum(X.^2, 1);
	X_index = (Xnorm2 > variance_Thresh);
           
	Ynorm2 = sum(Y.^2, 1);
	Y_index = (Ynorm2 > variance_Thresh);
		
    XY_index=X_index|Y_index;	
	X = X(:, XY_index);    
	Y = Y(:, XY_index);
    save DX_all X;
    save DY_all Y;
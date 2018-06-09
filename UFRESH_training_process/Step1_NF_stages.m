clear;
stage = 4;
%% ############  CHECK ME  ####################
if stage == 1
    savefilename_X = sprintf('../Python/DX_and_DY/DX_all_NF.mat');
    savefilename_Y = sprintf('../Python/DX_and_DY/DY_all_NF.mat');
    directory_x = 'TrainingData/NF_0'; %Stage 1
else
    savefilename_X = sprintf('../Python/DX_and_DY/DX_all_NF%i.mat',stage);
    savefilename_Y = sprintf('../Python/DX_and_DY/DY_all_NF%i.mat',stage);
    directory_x = sprintf('TrainingData/NF_%i',stage-1); %Stage 2+
end

tic
addpath('ksvd');
addpath('ksvd/ksvdbox');
addpath('ksvd/ksvdbox/private_ccode');
addpath('ksvd/ompbox');
addpath('utils');
addpath('Step2_Kmeans_clustering');


pattern = '*.bmp';
directory_y = 'TrainingData/GT'; 

XpathCell = glob(directory_x, pattern );
Xcell = load_images( XpathCell );

YpathCell = glob(directory_y, pattern );
Ycell = load_images( YpathCell );
blocksize = [5, 5]; % the size of each image patch.
% trainnum = 100000;
trainnum = 20000;
% trainnum = 100000;
variance_Thresh = 0.02;
X = [];
Y = [];

[n1, n2] = size(Xcell);
for i = 1: n1
    fprintf('%i out of %i\n',i,n1)        
%     X_temp = rot90(Xcell{i}, j-1);         
%     Y_temp = rot90(Ycell{i}, j-1);
    X_temp = Xcell{i};         
    Y_temp = Ycell{i};

    p = ndims(Y_temp); % Number of dimensions.
    ids = cell(p,1); % indices
    [ids{:}] = reggrid(size(Y_temp)-blocksize+1, trainnum, 'eqdist');
    Y_a = sampgrid(Y_temp,blocksize,ids{:});
    X_a = sampgrid(X_temp,blocksize,ids{:});

    X = [X X_a];            
    Y = [Y Y_a];
end
% Keep a copy of the non zero mean set
X1 = X;
Y1 = Y;
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
X1 = X1(:,XY_index);
Y1 = Y1(:,XY_index);
save(savefilename_X, 'X','X1');
save(savefilename_Y, 'Y','Y1');
toc
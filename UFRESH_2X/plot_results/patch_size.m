%% Average 5 + 14
clear
sz = [3:8,10];
psnrs = [32.279,32.446,32.521,32.521,32.519,32.500,32.464];
ssims = [0.9098,0.9126,0.9140,0.9141,0.9143,0.9142,0.9142];
traintime = [151,307,471,716,967,1218,1825]/60;
tpp = [12.3,15.0,22.5,22.7,27.5,38.3,54.3];
%% Set 5
clear
sz = [3:8,10];
psnrs = [35.120,35.319,35.429,35.433,35.422,35.411,35.366];
ssims = [0.9437,0.9462,0.9477,0.9479,0.9479,0.9480,0.9480];
traintime = [151,307,471,716,967,1218,1825]/60;
tpp = [12.3,15.0,22.5,22.7,27.5,38.3,54.3];
%% Set 14
clear
sz = [3:8,10];
psnrs = [31.265,31.420,31.482,31.481,31.482,31.461,31.428];
ssims = [0.8977,0.9006,0.9020,0.9021,0.9022,0.9022,0.9021];
traintime = [151,307,471,716,967,1218,1825]/60;
tpp = [12.3,15.0,22.5,22.7,27.5,38.3,54.3];

%% Plot psnr and ssim
for i = 1:length(sz)
    xticks{i} = sprintf('%ix%i',sz(i),sz(i));
end

figure;
plotme(sz,psnrs,'b',1);
title('PSNR');
xlabel('Patch Size')
ylabel('PSNR(dB)');
legend('230k Training Samples');
set_gca;
set(gca,'xtick',[3:8,10],'xticklabel',xticks);

figure;
plotme(sz,ssims,'b',1);
title('SSIM');
xlabel('Patch Size')
ylabel('SSIM');
legend('230k Training Samples');
set_gca;
set(gca,'xtick',[3:8,10],'xticklabel',xticks);


%% Plot train and test time
figure;
plotme(sz,traintime,'b',1);
title('Training Time for Different Patch Sizes');
xlabel('Patch Size')
ylabel('Training Time(mins)');
legend('230k Training Samples');
set_gca;
set(gca,'xtick',[3:8,10],'xticklabel',xticks);

figure;
plotme(sz,tpp,'b',1);
title('Encoding Time for Different Patch Sizes');
xlabel('Patch Size')
ylabel('Encoding Time(\mus/pixel)');
legend('8192 Centroids');
set_gca;
set(gca,'xtick',[3:8,10],'xticklabel',xticks);

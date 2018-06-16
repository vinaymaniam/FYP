clear
%% 1.17M samples
nctrds = [128,256,512,1024,2048,4096,8192,16384];
tcluster = [11.8,26.4,24.5,45.1,27.1,86.8,212,1221]/60;
tmap = [40,80,173,358,760,1446,2957,8112]/60;
trun = [13.01,13.42,14.49,16.12,16.77,18.83,23.79,28.79];

%% 230k samples
nctrds = [128,256,512,1024,2048,4096,8192,16384];
tcluster = [3.3,4.3,9.3,3.7,8,42.8,137,444.2]/60;
tmap = [7,14,28,60,119,239,471,916]/60;
trun = [13.01,13.42,14.49,16.12,16.77,18.83,23.79,28.79];

%% Plot timings
ttrain = tcluster + tmap;

% figure;
% plotme(nctrds,tmap,'r',1); hold on;
% plotme(nctrds,tcluster,'b',1); hold on;
% plotme(nctrds,ttrain,'k',1); hold on;
% title('Training Time for Different Numbers of Centroids')
% xlabel('Number of Centroids')
% ylabel('Time(minutes)')
% legend('Mapping Calculation','Clustering Time','Total Training Time')
% set_gca

figure;
plotme(nctrds,trun,'r',1); hold on;
plot(nctrds,11+0.15*sqrt(nctrds),'b:','LineWidth',2)
title('Super Resolution Time for Different Numbers of Centroids')
xlabel('Number of Centroids')
ylabel('Time per Pixel(\mus/pixel)')
legend('Encoding Time','0.15\surd{n}+15')
set_gca

%% Plot PSNRS 
psnr5_1M = [35.28,35.33,35.38,35.40,35.42,35.43,35.44,35.46];
psnr14_1M = [31.39,31.43,31.44,31.48,31.49,31.48,31.50,31.51];
psnr5 = [35.31,35.33,35.37,35.37,35.40,35.42,35.43,35.42];
psnr14 = [31.39,31.40,31.43,31.44,31.46,31.46,31.48,31.48];

figure;
plotme(nctrds,psnr5,'r',2); hold on;
plotme(nctrds,psnr5_1M,'b',2); hold on;
title('[Set 5] PSNR for Different Numbers of Centroids')
xlabel('Number of Centroids');
ylabel('PSNR(dB)')
legend('230k Training Samples','1.17M Training Samples');
set_gca

figure;
plotme(nctrds,psnr14,'r',2); hold on;
plotme(nctrds,psnr14_1M,'b',2); hold on;
title('[Set 14] PSNR for Different Numbers of Centroids')
xlabel('Number of Centroids');
ylabel('PSNR(dB)')
legend('230k Training Samples','1.17M Training Samples');
set_gca

%% Plot SSIM 
ssim5_1M = [0.9463,0.9464,0.9468,0.9473,0.9476,0.9475,0.9473,0.9473];
ssim14_1M = [0.9005,0.9007,0.9010,0.9016,0.9019,0.9018,0.9017,0.9018];
ssim5 = [0.9469,0.9467,0.9472,0.9474,0.9474,0.9475,0.9477,0.9475];
ssim14 = [0.9009,0.9009,0.9014,0.9016,0.9017,0.9018,0.9020,0.9018];

figure;
plotme(nctrds,ssim5,'r',2); hold on;
plotme(nctrds,ssim5_1M,'b',2); hold on;
title('[Set 5] SSIM for Different Numbers of Centroids')
xlabel('Number of Centroids');
ylabel('SSIM')
legend('230k Training Samples','1.17M Training Samples');
set_gca

figure;
plotme(nctrds,ssim14,'r',2); hold on;
plotme(nctrds,ssim14_1M,'b',2); hold on;
title('[Set 14] SSIM for Different Numbers of Centroids')
xlabel('Number of Centroids');
ylabel('SSIM')
legend('230k Training Samples','1.17M Training Samples');
set_gca






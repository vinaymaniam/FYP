clear
%% Plot PSNRS 
psnrs = [32.505,32.519,32.522,32.521,32.498,32.475,32.427];
ssims = [0.9137,0.9137,0.9138,0.9140,0.9138,0.9136,0.9134];
nn = [12,48,96,192,384,768,3072];

figure;
plotme(nn,psnrs,'b',2); hold on;
title('Average PSNR across Set 5 and 14')
xlabel('Neighbourhood Size');
ylabel('PSNR(dB)')
legend('230k Training Samples');
set_gca

figure;
plotme(nn,ssims,'b',2); hold on;
title('Average SSIM across Set 5 and 14')
xlabel('Neighbourhood Size');
ylabel('SSIM')
axis([-inf,inf,0.913,0.915])
legend('230k Training Samples');
set_gca






clear
%% Plot PSNRS 
sz = 1:8;

psnr_nf = [32.08668537,32.16170166,32.18397369,32.2043163,32.21539057,32.22443354,32.22544923,32.22990377];
psnr_f = [32.44911408,32.49526622,32.51344023,32.52072381,32.52669423,32.52854944,32.52903512,32.53058971];

ssim_nf = [0.907828599,0.909230396,0.909638344,0.909913306,0.910089949,0.910244456,0.910296226,0.910385508];
ssim_f = [0.912646971,0.91354108,0.913870848,0.914038565,0.914124191,0.91418578,0.914219798,0.914253046];

figure;
plotme(sz,psnr_nf,'b',1); hold on;
plotme(sz,psnr_f,'r',1); hold on;
rl1 = refline(0,max(psnr_nf));
rl1.Color = 'b';
rl1.LineStyle = '--';
rl2 = refline(0,max(psnr_f));
rl2.Color = 'r';
rl2.LineStyle = '--';
title('Average PSNR across Set 5 and 14')
xlabel('Ensemble Size');
ylabel('PSNR(dB)')
legend('Bicubic Input','FRESH Input');
set_gca


figure;
plotme(sz,ssim_nf,'b',1); hold on;
plotme(sz,ssim_f,'r',1); hold on;
rl1 = refline(0,max(ssim_nf));
rl1.Color = 'b';
rl1.LineStyle = '--';
rl2 = refline(0,max(ssim_f));
rl2.Color = 'r';
rl2.LineStyle = '--';
title('Average SSIM across Set 5 and 14')
xlabel('Ensemble Size');
ylabel('SSIM')
legend('Bicubic Input','FRESH Input');
set_gca





clear
varth = [0.003,0.01,0.02,0.025,0.03,0.04,0.06,0.1];

psnrs = [32.499,32.516,32.523,32.523,32.520,32.516,32.507,32.500];
ssims = [0.9141,0.9143,0.9143,0.9143,0.9143,0.9143,0.9143,0.9142];



figure;
plotme(varth,psnrs,'b',1);
title('Average PSNR Across Set 5 and 14');
xlabel('Variance Threshold')
ylabel('PSNR(dB)');
legend('Config 2.1');
set_gca;

figure;
plotme(varth,ssims,'b',1);
title('Average SSIM Across Set 5 and 14');
xlabel('Variance Threshold')
ylabel('SSIM');
axis([-inf,inf,-inf,max(ssims)+0.0001])
legend('Config 2.1');
set_gca;

clear
sz = [3:8,10];

psnrs = [32.279,32.446,32.521,32.521,32.519,32.500,32.464];
ssims = [0.9098,0.9126,0.9140,0.9141,0.9143,0.9142,0.9142];
traintime = [151,307,471,716,967,1218,1825]/60;
tpp = [12.3,15.0,22.5,22.7,27.5,38.3,54.3];


for i = 1:length(sz)
    xticks{i} = sprintf('%ix%i',sz(i),sz(i));
end

figure;
plotme(sz,psnrs,'b',1);
title('Average PSNR Across Set 5 and 14');
xlabel('Patch Size')
ylabel('PSNR(dB)');
legend('230k Training Samples');
set_gca;
set(gca,'xtick',[3:8,10],'xticklabel',xticks);

figure;
plotme(sz,ssims,'b',1);
title('Average SSIM Across Set 5 and 14');
xlabel('Patch Size')
ylabel('SSIM');
legend('230k Training Samples');
set_gca;
set(gca,'xtick',[3:8,10],'xticklabel',xticks);

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

function [outputArg1,outputArg2] = visualiseImprovement(LR,SR,HR)
    % Display the change between image 1 and 2(FRESH and UFRESH)
    % Block size 2*n x 2*n
    n = 64;

    figure(1)
    imshow(LR)   

    figure(1);
    [x,y] = ginput(1);
    idx1 = int32([max(1,y-n):min(size(LR,1),y+n-1)]);
    idx2 = int32([max(1,x-n):min(size(LR,1),x+n-1)]);
    figure(2);
    subplot(1,3,1);
    imshow(LR(idx1,idx2))
    title(sprintf('FRESH(PSNR %.2fdB)',psnr(LR(idx1,idx2),HR(idx1,idx2))));
    subplot(1,3,2);
    imshow(SR(idx1,idx2))
    title(sprintf('UFRESH(PSNR %.2fdB)',psnr(SR(idx1,idx2),HR(idx1,idx2))));
    subplot(1,3,3);
    imshow(HR(idx1,idx2))
    title('GT');
    truesize(size(LR(idx1,idx2))*7)
    waitforbuttonpress;   
end


function fresh_gui

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  NOTICE                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          GUI software for FRESH Single-Image Super-Resolution
% 
% The software implements the super-resolution algorithm described in:
% [1] X. Wei and P. L. Dragotti, "Sampling piecewise smooth signals and its
%     application to image up-sampling", in ICIP, Sep. 2015, pp. 4293?4297. 
% [2] X. Wei and P. L. Dragotti, "FRESH - FRI-based single-image 
%     super-resolution algorithm", in IEEE Transactions on Image Processing, 
%     vol. 25, no. 8, pp. 3723-3735, Aug. 2016.
% 
% Authors:      Xiaoyao Wei           ivy.wei@imperial.ac.uk
%               Matteo Maggioni       m.maggioni@imperial.ac.uk
%               Pier Luigi Dragotti   p.dragotti@imperial.ac.uk
% 
% Copyright (c) 2015-2016 Imperial College London.
% All rights reserved.
% 
% Disclaimer:
% This work should be used for nonprofit purposes only. Any unauthorized 
% use of these routines for industrial or profit-oriented activities is 
% expressively prohibited. By downloading and/or using any of these files, 
% you implicitly agree to all the terms of the limited license included in 
% this package (Legal_Notice.txt).

dwtmode('per','nodisp');
% Center=load ('Center2048.mat');
% Map=load ('Map2048.mat');
% Center=Center.Center;
% Map=Map.Map;
block_size = 64;
level = 2;
diag_reg = 1;
lin_map = 1;
high_res_in = 0;
profile = 0;
truesize_img = 1;

filename = '';
img = [];
hLine = [];
fLR = [];
fBic = [];
fFRESH = [];

if isunix
    pnlHeight = 64;
    offx = 0;
    offy = 0;
else
    pnlHeight = 72;
    offx = 15;
    offy = 37;
end

last_ups_block = [];

fGUI = figure('Visible','off','Units','Normalized',...
    'MenuBar','none','ToolBar','none',...
    'Name','FRESH - Single-Image Super-Resolution',...
    'NumberTitle','off',...
    'OuterPosition',[0 0 1 1]/1.5,...
    'SizeChangedFcn',@resize_fcn,...
    'CloseRequestFcn',@close_fcn);
movegui(fGUI,'center')

hPnl1 = uipanel('Parent', fGUI, 'FontSize',14, 'Units', 'Normalized');
hPnl2 = uipanel('Parent', fGUI, 'FontSize',14, 'Units', 'Normalized');
hPnl3 = uipanel('Parent', fGUI, 'FontSize',14, 'Units', 'Normalized', ...
    'Visible','off', 'Position', [.7 .52 .275 .45]);
hPnl4 = uipanel('Parent', fGUI, 'FontSize',14, 'Units', 'Normalized', ...
    'Visible','off', 'Position', [.7 .05 .275 .45]);
hAx1 = axes('Parent', hPnl1, 'Units', 'Normalized', ...
    'Position', [.025 .025 .95 .95]);

hPnl2.Title = 'Parameters';
uicontrol('Parent',hPnl2,'Style', 'text','Units','pixels',...
    'String', 'Block Size', 'HorizontalAlignment','left',...
    'Position', [10 20 100 25], ...
    'Tooltip','Size (in px) of the block to upsample.');
uicontrol('Parent',hPnl2,'Style', 'popup','Units','pixels',...
    'String', {'32x32','64x64','128x128'},...
    'Value',2, 'Position', [5 5 100 25],...
    'Callback', @set_block_size_fcn);

uicontrol('Parent',hPnl2,'Style', 'text','Units','pixels',...
    'HorizontalAlignment','left',...
    'String', 'Upsampling Factor',...
    'Position', [115 20 100 25]);
uicontrol('Parent',hPnl2,'Style', 'popup','Units','pixels',...
    'String', {'2X','4X','8X'},...
    'Value',1,...
    'Position', [110 5 100 25],...
    'Callback', @set_level_fcn);

uicontrol('Parent',hPnl2,'Style', 'checkbox',...
    'String','Diagonal Reg.',...
    'Value',1, 'Position', [230 24 100 22],...
    'Callback', @set_diag_reg_fcn, ...
    'Tooltip','Enable FRESH diagonal regularization.');

uicontrol('Parent',hPnl2,'Style', 'checkbox',...
    'String','Linear Mapping',...
    'Value',1, 'Position', [230 5 100 22],...
    'Callback', @set_lin_map_fcn, ...
    'Tooltip','Enable FRESH error correction by linear mapping.');

uicontrol('Parent',hPnl2,'Style', 'checkbox',...
    'String','Test Mode',...
    'Position', [340 24 100 22],...
    'Callback', @set_high_res_in_fcn, ...
    'Tooltip',['<html>Downsample input block such that the upsampled<br /> ',...
    'result will have the same size as the original input.</html>']);

uicontrol('Parent',hPnl2,'Style', 'checkbox',...
    'String','Fast Profile',...
    'Position', [340 5 100 22],...
    'Callback', @set_profile_fcn, ...
    'Tooltip','Speed-up FRESH algorithm.');

uicontrol('Parent',hPnl2,'Style', 'checkbox',...
    'String','Truesize',...
    'Value',1, 'Position', [430 5 100 22],...
    'Callback', @set_truesize_fcn, ...
    'Tooltip','Show processed blocks in truesize.');

hBtnNewImg = uicontrol('Parent',hPnl2,'Style', 'pushbutton',...
    'String','Load New Image',...
    'Position', [470 6 100 20],...
    'Callback', @load_new_image_fcn);

hBtnUpsImg = uicontrol('Parent',hPnl2,'Style', 'pushbutton',...
    'String','Upsample Image',...
    'Position', [470 29 100 20],...
    'Callback', @upsample_image_fcn);

hPnl3.Title = 'Input Block';
hAx2 = axes('Parent', hPnl3, 'Units', 'Normalized', ...
    'Position', [.025 .025 .95 .95]);
imshow(ones(block_size)*fGUI.Color(1),'Parent',hAx2)

hPnl4.Title = 'Output Block';
hAx3 = axes('Parent', hPnl4, 'Units', 'Normalized', ...
    'Position', [.025 .025 .95 .95]);
imshow(ones(block_size)*fGUI.Color(1),'Parent',hAx3)

load_new_image_fcn;

fGUI.Visible = 'on';

set(fGUI,'WindowButtonDownFcn',@upsample_block_fcn);
set(fGUI,'WindowButtonMotionFcn',@upsample_block_fcn);

function resize_fcn(~, ~)
    f_px = getpixelposition(fGUI);
    pnl2_prcn = pnlHeight/f_px(4);
    pnl1_prcn = 1-pnl2_prcn;
    hPnl1.Position = [.025 pnl2_prcn+0.07 .65 pnl1_prcn-0.1];
    hPnl2.Position = [.025 .05 .65 pnl2_prcn];
    buttons = {hBtnNewImg, hBtnUpsImg};
    posPnl2 = getpixelposition(hPnl2);
    for i=1:length(buttons)
        pos = buttons{i}.Position;
        pos(1) = posPnl2(3)-110;
        buttons{i}.Position = pos;
    end
end
function close_fcn(~, ~)
    delete(fLR);
    delete(fBic);
    delete(fFRESH);
    delete(fGUI);
end
function set_block_size_fcn(source, ~)
    value = source.String{source.Value};
    block_size = str2double(value(1:find(value=='x',1,'first')-1));
end
function set_level_fcn(source, ~) 
    level = str2double(source.String{source.Value}(1));
end
function set_diag_reg_fcn(source, ~)
    diag_reg = source.Value;
end
function set_lin_map_fcn(source, ~)
    lin_map = source.Value;
end
function set_high_res_in_fcn(source, ~)
    high_res_in = source.Value;
end
function set_profile_fcn(source, ~)
    profile = source.Value;
end
function set_truesize_fcn(source, ~)
    truesize_img = source.Value;
    if ~truesize_img
        delete(fLR);
        delete(fBic);
        delete(fFRESH);
    end
end
function load_new_image_fcn(~, ~)
    [filename,folder] = uigetfile(...
        {'*.jpg; *.jpeg; *.bmp; *.png; *.tif; *.tiff',...
        'Supported Files (*.jpg,*.bmp,*.png,*tif)'},...
        'Select the input image file');
    if ~ischar(filename) && isempty(img)
        error('No Image Selected.');
    end
    if ischar(filename)
        info = imfinfo([folder,filename]);
        bitdepth = info.BitDepth;
        if strcmpi(info.ColorType,'truecolor')
            bitdepth = bitdepth/3;
        end
        img = double(imread([folder,filename]))/(2^bitdepth-1)*255;
        imshow(uint8(img), 'Parent',hAx1)
        hLine = line('Visible','off','Color','y','LineWidth',2, 'Parent',hAx1);
        hPnl1.Title = sprintf('Input image (%s) - %dx%dpx - Click anywhere on the image', ...
            filename, size(img,1),size(img,2));
    end
end
function upsample_image_fcn(~, ~)
    N = block_size;
    block_size = size(img,1);
    src.CurrentAxes.CurrentPoint = [size(img,2) size(img,1)];
    evnt.EventName = 'WindowMousePress';
    upsample_block_fcn(src, evnt);
    block_size = N;
end
%deng
function [ I_rec ] = reconst(I)
blocksize=[5,5];
[HI,WI]=size(I); 
p=blocksize(1);
C=zeros(p*p,(HI-p+1),(WI-p+1));
I_rec=zeros(size(I));
count=I_rec;
for i=1:1:HI-p+1
    for j=1:1:WI-p+1
    C(:,i,j)=im2col(I(i:i+p-1,j:j+p-1),[p p]);
    end
end
for M=1:HI-p+1   
    for N=1:WI-p+1
        C_patch=C(:,M,N)-mean(C(:,M,N));
        if sum(C_patch.^2, 1)>0.1 % threshold 
        MSE_rough=sqrt(sum((Center-repmat(C_patch,1,size(Center,2))).^2));
        mse=sort(MSE_rough(:));
        t=find(MSE_rough<=mse(1),1); % find the most minimum MSE       
        avgp=Map(:,t);
        avgp=reshape(avgp,p*p,p*p);
        % multiply by the projection                    
        patch_rec=avgp*C_patch+mean(C(:,M,N));
        else
        patch_rec=C(:,M,N);
        end
        patch_rec=reshape(patch_rec,[p p]);
        % assemble into the reconstructed image
        I_rec(M:M+p-1,N:N+p-1)=I_rec(M:M+p-1,N:N+p-1)+patch_rec;
        count(M:M+p-1,N:N+p-1)= count(M:M+p-1,N:N+p-1)+1; 
    end
end
I_rec=I_rec./count;
end


%deng


function upsample_block_fcn(source, event)
    point = floor(source.CurrentAxes.CurrentPoint(1,1:2));
    hLine.XData = [point(1)-block_size, point(1), point(1), ...
        point(1)-block_size, point(1)-block_size]+0.5;
    hLine.YData = [point(2)-block_size, point(2)-block_size, ...
        point(2), point(2), point(2)-block_size]+0.5;
    hLine.Visible = 'off';
    if all(point>=block_size) && all(point<=[size(img,2) size(img,1)])
        hPnl3.Visible = 'on';
        hPnl4.Visible = 'on';
        hLine.Color = 'y';
        hLine.Visible = 'on';
        if ~isequal(last_ups_block,point)
            block = double(img(point(2)-block_size+1:point(2),...
                point(1)-block_size+1:point(1),:));
            block_lo = block;
            if high_res_in
                for i=1:log2(level)
                    block_lo = dwt2(block_lo,'haar')/2;
                end
            end
            imshow(uint8(block_lo), [0 255], 'Parent',hAx2)
            hPnl3.Title = sprintf('Low-Res  - %dx%dpx', ...
                size(block_lo,1), size(block_lo,2));
            hPnl1.Title = sprintf('Input image (%s) - %dx%dpx - x=%d y=%d', ...
                filename, size(img,1),size(img,2), ...
                point(1)-block_size, point(2)-block_size);
            if strcmp(event.EventName,'WindowMousePress')
                hPnl4.Title = 'Processing...';
                pause(eps),drawnow
                tic;
                block_up = fresh(block, log2(level), ...
                    diag_reg, lin_map, high_res_in, profile);
                imshow(uint8(block_up), [0 255], 'Parent',hAx2) % deng
                block_up=block_up./256; %deng
                block_ups=reconst(block_up); % deng
                elapsed_time = toc;
                block_ups = circshift(block_ups,[0 0]);
                algo = sprintf('FRESH (%.2fs)', elapsed_time);
                last_ups_block = point;
                % show truesize blocks
                if truesize_img
                    if isempty(fLR) || ~ishandle(fLR)
                        fLR = figure('Visible','off',...
                            'MenuBar','none','ToolBar','none',...
                            'Name','Low-Res','NumberTitle','off');
                    end
                    hAxLR = axes('Parent',fLR);
                    fLR.Visible = 'off';
                    posGUI = getpixelposition(fGUI);
                    imshow(uint8(block), [0 255], 'Parent',hAxLR) %deng
                   
                    truesize(fLR)
                    posLR = getpixelposition(fLR);
                    fLR.Position = [posGUI(1)+posGUI(3)+offx, ...
                        posGUI(2)+posGUI(4)-posLR(4), ...
                        posLR(3), posLR(4)];
                    if isempty(fBic) || ~ishandle(fBic)
                        fBic = figure('Visible','off',...
                            'MenuBar','none','ToolBar','none',...
                            'Name','Bicubic','NumberTitle','off');
                    end
                    hAxBic = axes('Parent',fBic);
                    fBic.Visible = 'off';
                    block_ups_bic = imresize(block_lo,level,'cubic');
                    imshow(uint8(block_ups_bic), [0 255], 'Parent',hAxBic)
                    truesize(fBic)
                    posBic = getpixelposition(fBic);
                    fBic.Position = [posGUI(1)+posGUI(3)+offx, ...
                        posGUI(2)+posGUI(4)-posLR(4)-posBic(4)-offy, ...
                        posBic(3), posBic(4)];
                    if isempty(fFRESH) || ~ishandle(fFRESH)
                        fFRESH = figure('Visible','off',...
                            'MenuBar','none','ToolBar','none',...
                            'Name','FRESH','NumberTitle','off');
                    end
                    hAxFRESH = axes('Parent',fFRESH);
                    fFRESH.Visible = 'off';
                    %imshow(uint8(block_ups), [0 255], 'Parent',hAxFRESH)
                    imshow(block_ups, 'Parent',hAxFRESH) % deng
                    truesize(fFRESH)
                    posFRESH = getpixelposition(fFRESH);
                    fFRESH.Position = [posGUI(1)+posGUI(3)+offx, ...
                        posGUI(2)+posGUI(4)-posLR(4)-posBic(4)-posFRESH(4)-offy*2, ...
                        posFRESH(3), posFRESH(4)];

                    fLR.Visible = 'on';
                    fBic.Visible = 'on';
                    fFRESH.Visible = 'on';
                end
            else
                block_ups = imresize(block_lo,level,'cubic');
                algo = 'Bicubic Interpolation';
            end
            %imshow(uint8(block_ups), [0 255], 'Parent',hAx3)
            imshow(block_ups, 'Parent',hAx3) % deng
            hPnl4.Title = sprintf('%s - %dx%dpx', algo, ...
                size(block_ups,1), size(block_ups,2));
        end
    else
        hLine.Visible = 'on';
        hLine.Color = 'r';
    end
    figure(fGUI);
end
       
end

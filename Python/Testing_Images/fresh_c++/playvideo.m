% (c) 2011-2013 Noiseless Imaging Oy - Tampere, Finland
% CONFIDENTIAL

function playvideo(varargin)
% playvideo.m     requires playback_control.m

%%
pause(eps)
clk_init=round(sum(100*clock));rand('seed',clk_init); uid1=(10.^((1:8)-1))*round(8*rand(8,1)+1);

global playback_mode_new
global frames_per_second
global gamma_corr
global truesizeRatio
global normalize01
global save_frame
global save_all
global colormap_counter
global jjj
global jjjRelative
global jloop
global enable_crop
global zoom_
global shiftc
global imageVsHist
global forceScreenUpdate
global permuteData
global rotateData
global cellIndexUpdate
global cellIndex
global numChannelsPerInput
global toggleAxesHist
global toggleAxes
global showHelp
global size_3
global FigFontSize
global showText
global GodMode
global popUpText
global showPopUpText
global showPopUpTextTimeOn
global popUpTextToggle
global popUpTextStyle


if iscell(varargin(1))&&nargin==1
    varargin=varargin{1};
end
if ~iscell(varargin)
    varargin={varargin};
end
cellIndex=1;
varnames=cell(size(varargin));
for ii=1:numel(varargin)
    varargin{ii}=single(varargin{ii});
    varnames{ii}=inputname(min(ii,nargin));
    if isempty(varnames{ii})
        varnames{ii}='input#';
    end
end

% defaults
frames_per_second=15;
gamma_corr=1;
% gamma_corr=1/2^(6/8);

colormap_counter=1;    %% colormap_name='gray';
% colormap_counter=11;    %% colormap_name='hot';
% colormap_counter=13;   %% colormap_name='thermal';
nColorBarRange=15;
normalize01=1;
save_all=0;
enable_crop=0;
imageVsHist=1;
FigFontSize=15;
% showText=0;
GodMode=false;
popUpTextToggle=true;
% savefileformat='.jpg';
% savefileformat='.bmp';
% savefileformat='.tif';  % this is 8-bit tiff
savefileformat='.tif 16bit';    % this is 16-bit tiff



%% THERE ARE NO USEFUL OPTIONS BEYOND THIS POINT


screenSize = get(0, 'ScreenSize'); % screenWidth = screenSize(3); screenHeight = screenSize(4);
if strcmp(savefileformat,'.tif 16bit')
    savefileformatExt='.tif';
else
    savefileformatExt=savefileformat;
end
[size_1 size_2 size_3]=size(varargin{cellIndex});
for jj=1:numel(varargin)
    min_z(jj)=min(varargin{jj}(:));
    max_z(jj)=max(varargin{jj}(:));
end
min_zAll=min(min_z);
max_zAll=max(max_z);
dimS=[1 3 2;2 3 1;3 2 1];
permuteData=0;
rotateData=0;
toggleAxes=0;
textBelowHandle=[];
toggleAxesHist=2;
% truesizeRatio=1;
truesizeRatio=0;
current_time=now;
playback_mode='play';
playback_mode_new='play';
showText=1;
showHelp=0;
save_frame=0;
how_many_skip=0;
jloop=3;
meanHist=[];
varHist=[];
frame_time=now;
forceScreenUpdate=0;
frame_time_old=current_time;
actual_fps=frames_per_second;
frames_per_second_old=frames_per_second;
imageVsHist_old=imageVsHist;
gamma_corr_old=gamma_corr;
colormap_counter_old=colormap_counter;
normalize01_old=normalize01;
time_to_wait=1/frames_per_second;
jjj=0;
jjj_old=0;
dimensionS=[1 2 3];
playback_mode_old='pause';
shiftc=[0 0 0 0];
old_shiftc=shiftc;
shiftRate=1;
cellIndexUpdate=0;
popUpTextPersistenceSeconds=1.5;
numChannelsPerInput=ones(1,numel(varargin));



if numel(varargin)>1
    if numel(varargin)>3
        helpTextTemp=',...,';
    elseif numel(varargin)>2 % numel(varargin)==3
        helpTextTemp=',2,';
    else %  numel(varargin)==2
        helpTextTemp=',';
    end
    helpTextTemp=['0,1',helpTextTemp,num2str(numel(varargin))];
    helpTextTemp=[helpTextTemp,' cycle through different cell inputs'];
else
    helpTextTemp='';
end
helpText={'KEYBOARD HELP';'(case sensitive)';'';helpTextTemp;...
    'p playback/pause   P reverse   +/- framerate   </> frame stepping';...
    'g/G/b  \gamma-correction    n/N normalization    c/C colormap';...
    'l rewind to first frame   L jump to middle frame';...
    'z/Z/X/x zoom in/out     ARROWS pan     k crop';...
    'r/R rotate frame     o/O permute dimensions';...
    't/T truesize ratio';...
    '{/} decrease/increase color channels';...
    'v/V toggle pop-up text';...
    'u/U cross-section  and  live histogram of current frame';...
    'a/A toggle axes     f/F font size';...
    's save frame    S save all frames';...
    'h toggle help screen    q quit';...
    };


showPopUpText=~false;
showPopUpTextTimeOn=now;
popUpText=cell(1);
popUpText{1}=varnames{cellIndex};
while ~strcmp(playback_mode,'exit')
    if jloop<3
        if ~ishandle(figuHandle)
            playback_mode_new='exit';
        end
    end
    if popUpTextToggle&&showPopUpText
        forceScreenUpdate=true;
    end
    
    playback_mode=playback_mode_new;
    if strcmp(playback_mode,'pause')
        how_many_skip=0;
        pause(eps);
    end
    if strcmp(playback_mode,'play')
        jjj=jjj+1;
    end
    if strcmp(playback_mode,'reverse')
        jjj=jjj-1;
    end
    if strcmp(playback_mode,'stepForward')
        jjj=jjj+1;
        how_many_skip=0;
        playback_mode_new='pause';
    end
    if strcmp(playback_mode,'stepBack')
        jjj=jjj-1;
        how_many_skip=0;
        playback_mode_new='pause';
        %         playback_mode='pause';
    end
    if frames_per_second==0;  %% this should never happen (because it is taken care of this already in playback_control.m), but it is left here just in case...
        playback_mode_new='pause';
        frames_per_second=1;
    end
    
    if permuteData||rotateData
        if permuteData
            title('PERMUTING DIMENSIONS')
            pause(eps);
            if permuteData==1
                dims=dimS(1,:);
                dimS=dimS([2 3 1],:);
            end
            if permuteData==-1
                [~, dims]=sort(dimS(3,:));
                dimS=dimS([3 1 2],:);
            end
            dimensionS=dimensionS(dims);
            for jj=1:numel(varargin)
                varargin{jj}=permute(varargin{jj},dims);
            end
            jjj=1;
            jjj_old=1;
            permuteData=0;
        else % rotateData
            if rotateData==-1
                for jj=1:numel(varargin)
                    varargin{jj}=permute(varargin{jj}(:,end:-1:1,:),[2 1 3]);
                end
                dimensionS(2)=-dimensionS(2);
                dimensionS=dimensionS([2 1 3]);
            elseif rotateData==1
                for jj=1:numel(varargin)
                    varargin{jj}=permute(varargin{jj}(end:-1:1,:,:),[2 1 3]);
                end
                dimensionS(1)=-dimensionS(1);
                dimensionS=dimensionS([2 1 3]);
            end
            rotateData=0;
        end
        frame_time_old=now;  how_many_skip=0;
        [size_1 size_2 size_3]=size(varargin{cellIndex});
        jloop=2;
    end
    
    
    if cellIndexUpdate
        
        if cellIndexUpdate==-1
            cellIndex=cellIndex+1;
        else
            cellIndex=cellIndexUpdate;
        end
        cellIndex=mod(cellIndex-1,numel(varargin))+1;
        if any([size_1 size_2 size_3]~= [size(varargin{cellIndex},1),size(varargin{cellIndex},2),size(varargin{cellIndex},3)]);
            sizeRates=[size(varargin{cellIndex},1),size(varargin{cellIndex},2),size(varargin{cellIndex},3)]./[size_1 size_2 size_3];
            crop_rectm(1:2)=round((crop_rectm(1:2)-1).*sizeRates([2 1])+1);
            crop_rectM(1:2)=round((crop_rectM(1:2)).*sizeRates([2 1]));
            crop_rectm(3:4)=crop_rectM-crop_rectm(1:2)+1;
            %             midPoint=((crop_rectm(1:2)+crop_rectm(3:4)-1)/2).*([size(varargin{cellIndex},2) size(varargin{cellIndex},1)]./[size_2 size_1])+1;
            %         crop_rectm=[1 1 size_2 size_1];
            %         crop_rectM=[size_2 size_1];
            
            
            [size_1 size_2 size_3]=size(varargin{cellIndex});
            
        end
        cellIndexUpdate=0;
        min_z(cellIndex)=min(varargin{cellIndex}(:));
        max_z(cellIndex)=max(varargin{cellIndex}(:));
        
        popUpText=cell(1);
        popUpText{1}=varnames{cellIndex};
        showPopUpText=true;
        popUpTextStyle=2;
        showPopUpTextTimeOn=now;
    end
    
    
    jjjRelative=mod(jjj-1,size_3)+1;
    
    if jloop==3
        figuHandle=figure('KeyPressFcn',@playback_control);
        axesHandle=gca;
        [size_1 size_2 size_3]=size(varargin{cellIndex});
        jloop=2;
    end
    if jloop==2
        crop_rectm=[1 1 size_2 size_1];
        crop_rectM=[size_2 size_1];
        jloop=1;
    end
    
    if zoom_~=0    %% ZOOMING
        zstep=max(1,round([crop_rectM(1:2)-crop_rectm(1:2)+1]/8));
        if zoom_==1;
            crop_rectMt=crop_rectM-zstep;
            crop_rectt(1:2)=crop_rectm(1:2)+zstep;
            for jjk=1:2
                if crop_rectMt(jjk)>=crop_rectt(jjk)
                    crop_rectM(jjk)=crop_rectMt(jjk);
                    crop_rectm(jjk)=crop_rectt(jjk);
                end
            end
        elseif zoom_==-1;
            crop_rectMt=crop_rectM+zstep;
            crop_rectt(1:2)=crop_rectm(1:2)-zstep;
            crop_rectM(1)=min(crop_rectMt(1),size_2);
            crop_rectM(2)=min(crop_rectMt(2),size_1);
            crop_rectm(1:2)=max(1,crop_rectt(1:2));
        end
        crop_rectm(3:4)=crop_rectM(1:2)-crop_rectm(1:2)+1;
        %         time_to_wait=0;
        jloop=1;
    end
    if any(shiftc)    %% PANNING
        if all(sign(shiftc)==sign(old_shiftc))
            shiftRate=shiftRate*1.1;
            shiftc=round(shiftc*shiftRate);
        else
            shiftRate=1;
        end
        
        crop_rectmBkp=crop_rectm;
        crop_rectMBkp=crop_rectM;
        crop_rectm(1:2)=max([1 1],min([size_2 size_1],crop_rectm(1:2)+shiftc(1:2)));
        crop_rectM(1:2)=max([1 1],min([size_2 size_1],crop_rectM(1:2)+shiftc(3:4)));
        crop_rectm(3:4)=crop_rectM(1:2)-crop_rectm(1:2)+1;
        if ~zoom_  %% this is to avoid frame_to_show shrinking when panning beyond the image boundary
            if crop_rectm(3)<crop_rectmBkp(3)  %H
                crop_rectm(1)=crop_rectmBkp(1);
                crop_rectM(1)=crop_rectMBkp(1);
                crop_rectm(3)=crop_rectmBkp(3);
            end
            if crop_rectm(4)<crop_rectmBkp(4)  %V
                crop_rectm(2)=crop_rectmBkp(2);
                crop_rectM(2)=crop_rectMBkp(2);
                crop_rectm(4)=crop_rectmBkp(4);
            end
        end
        old_shiftc=shiftc;
        shiftc=0;
        jloop=1;
        time_to_wait=-1/frames_per_second;
    else
        shiftRate=shiftRate/1.1;
    end
    forceScreenUpdate_old=forceScreenUpdate;
    if gcf==figuHandle&&(jloop==1||(jjj~=jjj_old||(size_3==1&&~strcmp(playback_mode,'pause')))||frames_per_second_old~=frames_per_second||gamma_corr_old~=gamma_corr||normalize01_old~=normalize01||colormap_counter_old~=colormap_counter||zoom_~=0||imageVsHist~=imageVsHist_old||forceScreenUpdate);
        zoom_=0;
        if numChannelsPerInput(cellIndex)==1
            frame_to_show=varargin{cellIndex}(crop_rectm(2):crop_rectM(2),crop_rectm(1):crop_rectM(1),jjjRelative);
        else
            frame_to_show=zeros(1+crop_rectM(2)-crop_rectm(2),1+crop_rectM(1)-crop_rectm(1),3);
            for jChannel=1:numChannelsPerInput(cellIndex)
                frame_to_show(:,:,jChannel)=(1/255)*varargin{mod(jChannel-1+cellIndex-1,numel(varargin))+1}(max(1,min(size(varargin{mod(jChannel-1+cellIndex-1,numel(varargin))+1},1),(crop_rectm(2):crop_rectM(2)))),max(1,min(size(varargin{mod(jChannel-1+cellIndex-1,numel(varargin))+1},2),(crop_rectm(1):crop_rectM(1)))),jjjRelative);
            end
        end
        %             disp([crop_rectm(2),crop_rectM(2),crop_rectm(1),crop_rectM(1)]);
        if normalize01&&(imageVsHist==1)  % do not normalize or gamma-correct when displaying histogram
            if normalize01==3
                minNorm01=min(frame_to_show(:));
                maxNorm01=max(frame_to_show(:));
            elseif normalize01==2
                minNorm01=min_z(cellIndex);
                maxNorm01=max_z(cellIndex);
            elseif normalize01==1
                minNorm01=min_zAll;
                maxNorm01=max_zAll;
            end
            frame_to_show=(frame_to_show-minNorm01)/(maxNorm01-minNorm01);
        else
            minNorm01=0;
            maxNorm01=1;
        end
        ColorBarRange=linspace(0,1,nColorBarRange);
        if gamma_corr~=1&&(imageVsHist==1)   % do not normalize or gamma-correct when displaying histogram
            frame_to_show=sign(frame_to_show).*abs(frame_to_show).^gamma_corr;
            if normalize01
                if normalize01==3
                    minNorm01b=min(frame_to_show(:));
                    maxNorm01b=max(frame_to_show(:));
                    frame_to_show=(frame_to_show-minNorm01b)/(maxNorm01b-minNorm01b);
                    ColorBarRange=ColorBarRange*(maxNorm01b-minNorm01b)+minNorm01b;
                elseif normalize01==2
                    %                     minNorm01b=min_z(cellIndex);
                    %                     maxNorm01b=max_z(cellIndex);
                    %                     frame_to_show=(frame_to_show-minNorm01b)/(maxNorm01b-minNorm01b);
                    %                     ColorBarRange=ColorBarRange*(maxNorm01b-minNorm01b)+minNorm01b;
                elseif normalize01==1
                    %nothing to do, because max and min are already 0 and 1
                end
            end
            ColorBarRange=sign(ColorBarRange).*abs(ColorBarRange).^(1/gamma_corr);
        end
        ColorBarRange=ColorBarRange*(maxNorm01-minNorm01)+minNorm01;
        
        if imageVsHist==1
            
            %             keyboard
            %                 imagesc(frame_to_show,[0 1]);
            %             image(frame_to_show,'CDataMapping','scaled','Parent',axesHandle)
            image(max(0,min(1,frame_to_show)),'CDataMapping','scaled','Parent',axesHandle)
            %             image(adapthisteq(frame_to_show),'CDataMapping','scaled','Parent',axesHandle)
            set(axesHandle,'CLim',[0 1])
            if truesizeRatio>=0
                screenFrameRatio=min(screenSize(4:-1:3)./[size(frame_to_show,1) size(frame_to_show,2)]);
                if screenFrameRatio>=(1/0.7+eps)
                    if truesizeRatio==0
                        truesizeRatio=floor(0.7*screenFrameRatio);
                    else
                        truesizeRatio=min(truesizeRatio,floor(0.7*screenFrameRatio));
                    end
                else
                    if truesizeRatio==0
                        truesizeRatio=1/ceil(1.4/screenFrameRatio);;
                    else
                        truesizeRatio=min(truesizeRatio,1/ceil(1.4/screenFrameRatio));
                    end
                end
                truesize(figuHandle,round([size(frame_to_show,1) size(frame_to_show,2)]*truesizeRatio));
                truesizeRatio=-truesizeRatio;
            end
        elseif imageVsHist==2
            nBins=max(2,round((numel(frame_to_show)/100)));
            [histValues binCenters]=hist(frame_to_show(:),nBins);
            meanHist=mean(frame_to_show(:)); varHist=var(frame_to_show(:));
            hhh=bar(axesHandle,binCenters,histValues,1);
            hold on
            binsRange=linspace(min_z(cellIndex),max_z(cellIndex),400);
            plot(binsRange,numel(frame_to_show)*(binCenters(2)-binCenters(1))*1/sqrt(2*pi*varHist)*exp(-(binsRange-meanHist).^2/(2*varHist)),'r','linewidth',1.5);
            hold off
            set(axesHandle,'xlim',[min_z(cellIndex) max_z(cellIndex)]);
        elseif imageVsHist==3
            image(log(abs(fftshift(fft2(frame_to_show)))),'CDataMapping','scaled','Parent',axesHandle)
            %             set(axesHandle,'CLim',[0 1])
        else
            plot(axesHandle,frame_to_show(ceil(end/2),:))
            
            
            if normalize01==3
                set(axesHandle,'ylim',[min(frame_to_show(ceil(end/2),:)),max(frame_to_show(ceil(end/2),:))]);
            elseif normalize01==2
                set(axesHandle,'ylim',[min_z(cellIndex),max_z(cellIndex)]);
            elseif normalize01==1
                set(axesHandle,'ylim',[min_zAll,max_zAll]);
            else
                set(axesHandle,'ylim',[0,1]);
            end
            
        end
        imageVsHist_old=imageVsHist;
        
        if colormap_counter_old~=colormap_counter||jloop==1
            colormap_counter=mod(colormap_counter-1,14)+1;
            if colormap_counter<3;
                colormap_name='gray';
            elseif colormap_counter<5;
                colormap_name='bone';
            elseif colormap_counter<7;
                colormap_name='pink';
            elseif colormap_counter<9;
                colormap_name='jet';
            elseif colormap_counter<11;
                colormap_name='hsv';
            elseif colormap_counter<13;
                colormap_name='hot';
            elseif colormap_counter<15;
                colormap_name='custom_cMap';
                cMap = [ 0.0 0.0 0.0;  0.3 0.0 0.7;  1.0 0.2 0.0;  1.0 1.0 0.0;  1.0 1.0 1.0];
                cMapLabel='thermal';
            end
            colormap_counter_old=colormap_counter;
            if strcmp(colormap_name,'custom_cMap')
                % cMap and cMapLabel already defined above
            else
                cMap=colormap(colormap_name);
                cMapLabel=colormap_name;
            end
            if ~mod(colormap_counter,2);
                cMap=colormap(flipud(cMap));
                cMapLabel=['-',cMapLabel];
            end
            
            % makes colormap of 256 levels, instead of the usual 64
            NcMapExt=256;clear cMapExt;for jc=1:3;cMapExt(:,jc)=interp1(1:size(cMap,1),cMap(:,jc),linspace(1,size(cMap,1),NcMapExt));end
            cMap=cMapExt;
            colormap(cMap);
        end
        if jloop==1
            if enable_crop                 % Cropping
                truesize
                title('---- PLEASE CROP REGION FROM IMAGE  ----')
                hold on
                [dummy,crop_rectm]=imcrop(round(frame_to_show*255),cMap);
                hold off
                crop_rectm(1:2)=max(1,ceil(crop_rectm(1:2)));
                if crop_rectm(1)>size_2
                    crop_rectm(1)=max(1,crop_rectm(1)-size_2);
                end
                crop_rectm(3:4)=min([size_2 size_1]-crop_rectm(1:2)+1,round(crop_rectm(3:4)));
                crop_rectM=[crop_rectm(1)+crop_rectm(3)-1, crop_rectm(2)+crop_rectm(4)-1];
                time_to_wait=1/frames_per_second;
                frame_time_old=now;
                enable_crop=0;
            else
                jloop=0;
            end
            colormap(cMap)
        end
        if imageVsHist==1
            if toggleAxes
                axis on
            else
                axis off
            end
            if toggleAxes==2
                grid on
            else
                grid off
            end
            axis equal tight
            if toggleAxes
                set(axesHandle,'xticklabel',get(axesHandle,'xtick')+crop_rectm(1)-1);
                set(axesHandle,'yticklabel',get(axesHandle,'ytick')+crop_rectm(2)-1);
            end
        else
            if toggleAxesHist
                axis on
            else
                axis off
            end
            if toggleAxesHist==2
                grid on
            else
                grid off
            end
        end
        set(axesHandle,'fontsize',FigFontSize)
        if imageVsHist==1
            for jjn=1:nColorBarRange %% this is the colorbar!
                colorbarTextLabel=sprintf('%0.4g',ColorBarRange(jjn));
                colorbarTextLabel=[' ',colorbarTextLabel];
                textSideRightHandle=text(0.5+2*(crop_rectM(1)-crop_rectm(1)+1)/2,0.5+(1-(jjn-1+1)/(nColorBarRange-1+2))*(1+crop_rectM(2)-crop_rectm(2)),colorbarTextLabel,'HorizontalAlignment','left','VerticalAlignment','middle','color',cMap(1+round((size(cMap,1)-1)*(jjn-1)/(nColorBarRange-1)),:),'FontWeight','bold','fontsize',FigFontSize,'interpreter','none');
            end
            if showText
                textSideLeft={{['cMap=',num2str(colormap_counter)];cMapLabel};[varnames{cellIndex},' ',num2str(cellIndex),'/',num2str(numel(varargin))];['dims=[',num2str(dimensionS(1)),' ',num2str(dimensionS(2)),' ',num2str(dimensionS(3)),']'];['\gamma=',num2str(gamma_corr)];['norm01=',num2str(normalize01)];['truesizeRatio=',num2str(max(1,abs(truesizeRatio))),'/',num2str(max(1,abs(1/truesizeRatio)))]};
                for jj=1:numel(textSideLeft)
                    textSideLeftHandle=text(0.5,0.5+((jj)/(numel(textSideLeft)+1))*(1+crop_rectM(2)-crop_rectm(2)),textSideLeft{jj},'HorizontalAlignment','right','VerticalAlignment','middle','fontsize',FigFontSize,'interpreter','none');
                end
            end
            if showHelp
                for jj=1:numel(helpText)
                    helpTextHandle=text(0.5+0.5*(1+crop_rectM(1)-crop_rectm(1)),0.5+((jj)/(numel(helpText)+1))*(1+crop_rectM(2)-crop_rectm(2)),helpText{jj},'HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor',[0 0 0],'Color',[1 1 1],'FontSize',FigFontSize,'FontName','FixedWidth','interpreter','none');
                end
            end
            
            if popUpTextToggle&&showPopUpText
                popUpTextLast=popUpText;
                popUpTextPersistenceSecondsLeft=popUpTextPersistenceSeconds-(now-showPopUpTextTimeOn)*86400;
                if popUpTextPersistenceSecondsLeft>0
                    for jj=1:numel(popUpTextLast)
                        if popUpTextStyle==1
                            popUpTextSize=FigFontSize*4;%*sin(pi/2*popUpTextPersistenceSecondsLeft/popUpTextPersistenceSeconds);
                            shadowShift=0.007*mean(crop_rectm(3:4))*popUpTextSize/(FigFontSize*4);
                            helpTextHandle=text(shadowShift+0.5+0.5*(1+crop_rectM(1)-crop_rectm(1)),shadowShift+0.5+((jj)/(numel(popUpTextLast)+1))*(1+crop_rectM(2)-crop_rectm(2)),popUpTextLast{jj},'HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor','none','Color',[0 0 0],'FontSize',popUpTextSize,'FontName','FixedWidth','interpreter','none','margin',1,'FontWeight','bold');
                            helpTextHandle=text(0.5+0.5*(1+crop_rectM(1)-crop_rectm(1)),0.5+((jj)/(numel(popUpTextLast)+1))*(1+crop_rectM(2)-crop_rectm(2)),popUpTextLast{jj},'HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor','none','Color',[1 1 1],'FontSize',popUpTextSize,'FontName','FixedWidth','interpreter','none','margin',1,'FontWeight','bold');
                        elseif popUpTextStyle==2
                            popUpTextSize=FigFontSize*4*sin(pi/2*popUpTextPersistenceSecondsLeft/popUpTextPersistenceSeconds);
                            shadowShift=0.007*mean(crop_rectm(3:4))*popUpTextSize/(FigFontSize*4);
                            helpTextHandle=text(shadowShift+0.5+0.5*(1+crop_rectM(1)-crop_rectm(1)),shadowShift+0.5+((jj)/(numel(popUpTextLast)+1))*(1+crop_rectM(2)-crop_rectm(2)),popUpTextLast{jj},'HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor','none','Color',[0 0 0],'FontSize',popUpTextSize,'FontName','FixedWidth','interpreter','none','margin',1,'FontWeight','bold');
                            helpTextHandle=text(0.5+0.5*(1+crop_rectM(1)-crop_rectm(1)),0.5+((jj)/(numel(popUpTextLast)+1))*(1+crop_rectM(2)-crop_rectm(2)),popUpTextLast{jj},'HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor','none','Color',[1 1 1],'FontSize',popUpTextSize,'FontName','FixedWidth','interpreter','none','margin',1,'FontWeight','bold');
                        else
                            helpTextHandle=text(0.5+0.5*(1+crop_rectM(1)-crop_rectm(1)),0.5+((jj)/(numel(popUpTextLast)+1))*(1+crop_rectM(2)-crop_rectm(2)),popUpTextLast{jj},'HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor',[0 0 0],'Color',[1 1 1],'FontSize',FigFontSize*2,'FontName','FixedWidth','interpreter','none','margin',FigFontSize);
                        end
                    end
                else
                    showPopUpText=false;
                end
            end
            
        else
            pause(eps)
        end
        if ~save_all
            if how_many_skip
                skippingText=' x';
            else
                skippingText=' v';
            end
            textAboveString=[playback_mode,':  frame# ',num2str(jjjRelative),'/',num2str(size_3),'  FPS set=',num2str(frames_per_second),' (actual=',num2str(actual_fps,'%2.1f'),skippingText,')'];
            textBelowString=['press H for help, Q to quit'];
            if imageVsHist==1
                if showText
                    textAboveHandle=text(0.5+0.5*(1+crop_rectM(1)-crop_rectm(1)),0.5,textAboveString,'HorizontalAlignment','center','VerticalAlignment','bottom','fontsize',FigFontSize,'interpreter','none');
                    textBelowHandle=text(0.5+0.5*(1+crop_rectM(1)-crop_rectm(1)),0.5+(1+crop_rectM(2)-crop_rectm(2)),textBelowString,'HorizontalAlignment','center','VerticalAlignment','top','fontsize',FigFontSize,'interpreter','none');
                end
            elseif imageVsHist==2;
                if showText
                    title(textAboveString);
                end
                xlabel(['mean = ',num2str(meanHist),'  st.dev = ',num2str(sqrt(varHist))]);
            elseif imageVsHist==0;
                if showText
                    title(textAboveString);
                end
            end
        else
            textAboveHandle=text(0.5+0.5*(1+crop_rectM(1)-crop_rectm(1)),0.5,['   saving frame# ',num2str(jjjRelative),'/',num2str(size_3)],'HorizontalAlignment','center','VerticalAlignment','bottom','fontsize',FigFontSize,'interpreter','none');
            textBelowHandle=text(0.5+0.5*(1+crop_rectM(1)-crop_rectm(1)),0.5+(1+crop_rectM(2)-crop_rectm(2)),'*** saving ***','HorizontalAlignment','center','VerticalAlignment','top','fontsize',FigFontSize,'interpreter','none');
        end
        pause(eps)
        frame_time=now;
        actual_fps=(1+how_many_skip)/(frame_time-frame_time_old)/86400;
        if strcmp(playback_mode,'reverse')||strcmp(playback_mode,'play')   %% this tries to maintain (if at all possible) the desired frame rate
            if strcmp(playback_mode_old,'reverse')||strcmp(playback_mode_old,'play')
                time_to_wait=time_to_wait+0.5*(1/frames_per_second-1/actual_fps);
            else
                time_to_wait=1/frames_per_second;
                actual_fps=frames_per_second;
            end
            time_to_waitt=time_to_wait*~save_all;
            if (time_to_wait<(-1/frames_per_second))&&~(save_frame||save_all)
                how_many_skip=round(-time_to_wait*frames_per_second);
                if strcmp(playback_mode,'play')&&strcmp(playback_mode_new,'play')
                    jjj=jjj+how_many_skip;
                elseif strcmp(playback_mode,'reverse')&&strcmp(playback_mode_new,'reverse')
                    jjj=jjj-how_many_skip;
                end
                jjjRelative=mod(jjj-1,size_3)+1;
            else
                how_many_skip=0;
            end
            while (strcmp(playback_mode,'reverse')||strcmp(playback_mode,'play'))&&time_to_waitt>0
                pause(max(0,min(1/25,time_to_waitt)))
                time_to_waitt=time_to_waitt-1/25;
            end
        end
        frame_time_old=frame_time;
        frames_per_second_old=frames_per_second;
        gamma_corr_old=gamma_corr;
        %         playback_mode_old=playback_mode;
        normalize01_old=normalize01;
    end
    if strcmp(playback_mode,'reverse')||strcmp(playback_mode,'play')
    else
        pause(1/100)
    end
    jjj_old=jjj;
    if (save_frame||save_all)&&jjjRelative>0
        clk_init=round(sum(100*clock));rand('seed',clk_init);    uid2=(10.^((1:8)-1))*round(8*rand(8,1)+1);
        if save_all
            uid2=0;
        end
        frame_number_max_digits=ceil(log10(size_3+1));
        frame_number_digits=ceil(log10(jjjRelative+1));
        save_frame_filename=['frame_',repmat('0',[1 frame_number_max_digits-frame_number_digits]),num2str(jjjRelative),'_I',num2str(cellIndex),'D',num2str(dimensionS*[100 10 1]'),'___uid',num2str(uid1),num2str(uid2)];
        if ~exist([save_frame_filename,savefileformatExt],'file')
            if strcmp(savefileformat,'.tif 16bit')
                imwrite(uint16(frame_to_show*(2^16-1)),[save_frame_filename,savefileformatExt]);
            else
                imwrite(frame_to_show,[save_frame_filename,savefileformatExt]);
            end
            if ishandle(textBelowHandle),delete(textBelowHandle),end
            textBelowHandle=text('Interpreter','none','Position',[0.5+0.5*(1+crop_rectM(1)-crop_rectm(1)),2.5+crop_rectM(2)-(crop_rectm(2)+1)],'String',['frames saved as  ',save_frame_filename,savefileformat],'HorizontalAlignment','center','VerticalAlignment','top','fontsize',FigFontSize);
        else %% file already exist, thus assumes that saving is complete.
            if save_all
                save_all=0;
                %                 playback_mode='pause';
                playback_mode_new='pause';
                if ishandle(textAboveHandle),delete(textAboveHandle),end
                textAboveString=[playback_mode,':  frame# ',num2str(jjjRelative),'/',num2str(size_3),'  FPS set=',num2str(frames_per_second),' (actual=',num2str(actual_fps,'%2.1f'),skippingText,')'];
                textAboveHandle=text(0.5+0.5*(1+crop_rectM(1)-crop_rectm(1)),0.5,textAboveString,'HorizontalAlignment','center','VerticalAlignment','bottom','fontsize',FigFontSize,'interpreter','none');
                save_frame_filename=['frame_',repmat('X',[1 frame_number_max_digits]),'_I',num2str(cellIndex),'D',num2str(dimensionS*[100 10 1]'),'___uid',num2str(uid1),num2str(uid2)];
            end
        end
        if ishandle(textBelowHandle),delete(textBelowHandle),end
        if save_frame
            textBelowHandle=text('Interpreter','none','Position',[0.5+0.5*(1+crop_rectM(1)-crop_rectm(1)),2.5+crop_rectM(2)-(crop_rectm(2)+1)],'String',['frame saved as  ',save_frame_filename,savefileformat],'HorizontalAlignment','center','VerticalAlignment','top','FontWeight','bold','fontsize',FigFontSize,'interpreter','none');
        else
            textBelowHandle=text('Interpreter','none','Position',[0.5+0.5*(1+crop_rectM(1)-crop_rectm(1)),2.5+crop_rectM(2)-(crop_rectm(2)+1)],'String',['all frames saved as  ',save_frame_filename,savefileformat],'HorizontalAlignment','center','VerticalAlignment','top','FontWeight','bold','fontsize',FigFontSize,'interpreter','none');
        end
        save_frame=0;
    end
    if forceScreenUpdate_old
        forceScreenUpdate=0;
        forceScreenUpdate_old=0;
    end
    playback_mode_old=playback_mode;
    
    if strcmp(playback_mode,'stepBack')||strcmp(playback_mode,'stepForward')
        forceScreenUpdate=1;
    end
    if GodMode
        GodMode=false;
        keyboard
    end
end
if strcmp(playback_mode,'exit')
    if ishandle(figuHandle)
        set(figuHandle,'KeyPressFcn','')
        if ishandle(textBelowHandle),delete(textBelowHandle),end
    end
end
if 0
    if ishandle(figuHandle)
        close(figuHandle)
    end
end
pause(eps)

end

%% playback control

% (c) 2011-2012 Noiseless Imaging Oy - Tampere, Finland
% CONFIDENTIAL

function playback_control(src,evnt)

global playback_mode_new
global save_all
global frames_per_second
global gamma_corr
global normalize01
global save_frame
global colormap_counter
global showHelp
global jjj
global jjjRelative
global jloop
global truesizeRatio
global enable_crop
global zoom_
global shiftc
global imageVsHist
global forceScreenUpdate
global permuteData
global rotateData
global cellIndexUpdate
global toggleAxes
global toggleAxesHist
global FigFontSize
global size_3
global showText
global GodMode
global cellIndex
global numChannelsPerInput
global popUpText
global showPopUpText
global showPopUpTextTimeOn
global popUpTextToggle
global popUpTextStyle

% if ~isempty(evnt.Key)

if strcmp(evnt.Key,'uparrow')
    shiftc=[0 -1 0 -1];
end
if strcmp(evnt.Key,'downarrow')
    shiftc=[0 1 0 1];
end
if strcmp(evnt.Key, 'leftarrow')
    shiftc=-[1 0 1 0];
end
if strcmp(evnt.Key,'rightarrow')
    shiftc=[1 0 1 0];
end

if isempty(evnt.Character), return; end


if evnt.Character == 'p'
    if strcmp(playback_mode_new,'pause')
        playback_mode_new='play';
        frames_per_second=max(1,frames_per_second);
    else
        playback_mode_new='pause';
        forceScreenUpdate=1;
    end
end
if evnt.Character == 'P'
    if strcmp(playback_mode_new,'reverse')
        playback_mode_new='play';
    else
        playback_mode_new='reverse';
    end
    frames_per_second=max(1,frames_per_second);
end

if ~isempty(str2num(evnt.Character))
    if str2num(evnt.Character)==round(abs(str2num(evnt.Character))); %% avoid strange things, e.g., imaginary unit when pressing 'I' key
        cellIndexUpdate=str2num(evnt.Character);
        if cellIndexUpdate==0;
            cellIndexUpdate=-1;
        end
        forceScreenUpdate=1;
    end
end

if evnt.Character == '>'||evnt.Character == '.'
    playback_mode_new='stepForward';
end
if evnt.Character == '<'||evnt.Character == ','
    playback_mode_new='stepBack';
end

if evnt.Character == 'l'
    jjj=0;
    playback_mode_new='stepForward';
end

if evnt.Character == 'L'
    jjj=floor(size_3/2);
    playback_mode_new='stepForward';
end

if evnt.Character == '#'
    GodMode=true;
end



if evnt.Character == '+'||evnt.Character == '='
    frames_per_second=frames_per_second+1;
end
if evnt.Character == '-'||evnt.Character == '_'
    frames_per_second=frames_per_second-1;
    frames_per_second=max(0,frames_per_second);
end

if evnt.Character == 't'
    truesizeRatio=max(eps,abs(truesizeRatio));
    if truesizeRatio>=1
        truesizeRatio=round(truesizeRatio+1);
    else
        truesizeRatio=1/max(1,round((1/truesizeRatio)-1));
    end
    forceScreenUpdate=1;
end
if evnt.Character == 'T'
    truesizeRatio=max(eps,abs(truesizeRatio));
    if truesizeRatio>=2
        truesizeRatio=round(truesizeRatio-1);
    else
        truesizeRatio=1/max(1,round((1/truesizeRatio)+1));
    end
    forceScreenUpdate=1;
end


if frames_per_second==0;
    playback_mode_new='pause';
    forceScreenUpdate=1;
    frames_per_second=1;
end


if evnt.Character == 'f'
    FigFontSize=FigFontSize-1;
    FigFontSize=max(1,FigFontSize);
    forceScreenUpdate=1;
end
if evnt.Character == 'F'
    FigFontSize=FigFontSize+1;
    forceScreenUpdate=1;
end


if evnt.Character == 'g'
    gamma_corr=gamma_corr/2^(1/8);
end
if evnt.Character == 'G'
    gamma_corr=gamma_corr*2^(1/8);
    gamma_corr=max(eps,gamma_corr);
end
if evnt.Character == 'b'||evnt.Character == 'B'
    gamma_corr=1;
end


if evnt.Character == 'h'||evnt.Character == 'H'
    showHelp=xor(showHelp,imageVsHist==1);
    imageVsHist=1;
    forceScreenUpdate=1;
end


if evnt.Character == 'n'
    normalize01=mod(normalize01+1,4);
elseif evnt.Character == 'N'
    normalize01=mod(normalize01-1,4);
end

if evnt.Character == 'k'||evnt.Character == 'K'
    jloop=2;
    enable_crop=1;
end




if evnt.Character == 'c'
    colormap_counter=colormap_counter+1;
end
if evnt.Character == 'C'
    colormap_counter=colormap_counter-1;
end


if evnt.Character == 's'
    save_frame=1;
    playback_mode_new='pause';
    %     forceScreenUpdate=1;  %% this is not enabled, in order to keep visible the text about filename
    imageVsHist=1;
end
if evnt.Character == 'S'
    save_all=1;
    jjj=0;
    playback_mode_new='play';
    imageVsHist=1;
end


if evnt.Character == 'a'
    if imageVsHist==1
        toggleAxes=mod(toggleAxes+1,3);
    else
        toggleAxesHist=mod(toggleAxesHist+1,3);
    end
    forceScreenUpdate=1;
end
if evnt.Character == 'A'
    if imageVsHist==1
        toggleAxes=mod(toggleAxes-1,3);
    else
        toggleAxesHist=mod(toggleAxesHist-1,3);
    end
    forceScreenUpdate=1;
end


if evnt.Character == 'x'||evnt.Character == 'Z'
    zoom_=-1;
end
if evnt.Character == 'z'||evnt.Character == 'X'
    zoom_=1;
end


if evnt.Character =='q'||evnt.Character == 'Q'
    playback_mode_new='exit';
end

if evnt.Character =='e'||evnt.Character == 'E'
    showText=~showText;
    forceScreenUpdate=1;
end

if evnt.Character =='v'||evnt.Character == 'V'
    popUpTextToggle=~popUpTextToggle;
end



if evnt.Character =='U'
    imageVsHist=mod(imageVsHist+1,4);
    forceScreenUpdate=1;
end
if evnt.Character == 'u'
    imageVsHist=mod(imageVsHist-1,4);
    forceScreenUpdate=1;
end

if evnt.Character == 'o'
    permuteData=1;
    forceScreenUpdate=1;
end
if evnt.Character =='O'
    permuteData=-1;
    forceScreenUpdate=1;
end

if evnt.Character == 'r'
    rotateData=1;
    forceScreenUpdate=1;
end
if evnt.Character =='R'
    rotateData=-1;
    forceScreenUpdate=1;
end

if evnt.Character =='{'
    numChannelsPerInput(cellIndex)=max(1,min(3,numChannelsPerInput(cellIndex)-1));
    forceScreenUpdate=1;
end
if evnt.Character =='}'
    numChannelsPerInput(cellIndex)=max(1,min(3,numChannelsPerInput(cellIndex)+1));
    forceScreenUpdate=1;
end

showPopUpText=true;
popUpText=cell(1);
popUpText{1}=evnt.Character;
popUpTextStyle=1;
showPopUpTextTimeOn=now;

jjjRelative=mod(jjj-1,size_3)+1;
%
% if strcmp(evnt.Key,'F1')  %% does not work!
%     disp('f1');
% end
%
% disp(evnt.Key)



end

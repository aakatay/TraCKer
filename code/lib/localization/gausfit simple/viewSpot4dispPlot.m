
% data
load spotInfo;
numFrmDisp = 1000;
numFrmDisp = nSpots;
cropDisp;  % keeps the spots to be displayed

% to be removed
err2by2 = err2G5by5_2*0+1;


y5y1=x5x1;y5y2=x5x2;
[n,m,m_]=size(Ifit2);
% define output XY array
[iyTrack ixTrack] = find(xBW>0); % positions of spot data in the X array

frstSpot = 1;
ixIn = frstSpot:SP;
spPrev = 0;

% inactive options
isGenSpot = 0;
intPeak = intPeak_;

numSpots = length(isFit5by5);

for i = 1 : length(errCenter)
    if isempty(cell2mat(errCenter(i)))
    errCenter_(i) = nan;
    else
    errCenter_(i) = 1;        
    end
end

%intFit(200)=0;
%err2by2plot(200)=0;

tit2 = 'controls';
hFigCont = figure('DoubleBuffer','on','Menubar','none','Name',tit2,'NumberTitle','off','Colormap',gray(256));
figPos = [1921 -461 1080 1844];
set(hFigCont,'Position', figPos);        

% error plots
nFg = 7;
err2G5by5_2(find(err2G5by5_2==0)) = nan;
hText_ = uicontrol('style','text','BackgroundColor',[1 1 1],'String',sprintf('max. intPeak:%.02f',max(intPeak(:))),'Position',[900 1730 160 20]);
mxIntPeak = max(intPeak(:));
intPeak = intPeak/mxIntPeak;


mn = min(err2by2(:)); mx = max(err2by2(:)); err2by2plot = (err2by2-mn)/(mx-mn);

hText_ = uicontrol('style','text','BackgroundColor',[1 1 1],'String',sprintf('err2by2: %.02f - %.02f',mn,mx),'Position',[100 1500 160 20]);        
mn = min(err2G5by5_2(:)); mx = max(err2G5by5_2(:)); err2G5by5_2plot = (err2G5by5_2-mn)/(mx-mn);
mnErr = mn; mxErr = mx; 

hText_ = uicontrol('style','text','BackgroundColor',[1 1 1],'String',sprintf('err2G5by5_2: %.02f - %.02f',mn,mx),'Position',[800 1500 160 20]);        

hScrollButton = uicontrol('style','pushbutton','String','start scroll','Position',[100 1780 100 20],'Callback', 'updScrollState(hScrollButton,hTextSpotNumState)');     
hTextSpotNumState = uicontrol('style','text','BackgroundColor',[1 1 1],'String',sprintf('select spot'),'Position',[250 1780 100 20]);     
hTextSpotNum = uicontrol('style','text','BackgroundColor',[1 1 1],'String',sprintf('spot#: %i',0),'Position',[400 1780 100 20]);     
hSpotSlider = uicontrol('style','slider','units','pixel','Position',[100 1750 900 20], 'SliderStep',[1/numSpots 10/numSpots]);
listenSpot = addlistener(hSpotSlider,'ActionEvent',@(hObject, event) updSpotNum(hObject, event,hTextSpotNum,numSpots));
%listenScrollBtn = addlistener(hScrollButton,'ActionEvent',@(hObject, event) updScrollState(hObject, event,hScrollButton,hTextSpotNumState));

% clear out non proc. spots
mrkSz = 2;
subplot(nFg,1,1);   s = ixIn; % PLOT #1
[hAxInt,hLine1,hLine2] = plotyy(s,intPeak(ixIn),s,intFit(ixIn)); hold on;
set(hLine1,'Color',[0 0 1],'MarkerEdgeColor','k','MarkerSize',mrkSz,'LineStyle','.','Marker','*')
set(hLine2,'Color',[0 0 0],'MarkerEdgeColor','r','MarkerSize',mrkSz,'LineStyle','.','Marker','*')
legend('intPeak','intFit');
ylim([0 1])
ylabel(hAxInt(1),'intPeak') % left y-axis
ylabel(hAxInt(2),'intFit') % right y-axis
title('select intensity threshold to filter spots')
xlabel('spots');grid minor

subplot(nFg,1,2); % PLOT #2
[hAxErr,hLine1,hLine2] = plotyy(s,err2G5by5_2plot(ixIn),s,err2by2plot(ixIn)); hold on;
%plot(s,err2G5by5_2);
set(hLine1,'Color',[0 0 1],'MarkerEdgeColor','b','MarkerSize',mrkSz,'LineStyle','.','Marker','*')
set(hLine2,'Color',[0 0 0],'MarkerEdgeColor','g','MarkerFaceColor',[0.5 0.5 1],'MarkerSize',mrkSz,'LineStyle','.','Marker','*')        
legend('err2G5by5\_2','err2by2');
ylabel(hAxErr(1),'err2G5by5\_2') % left y-axis
ylabel(hAxErr(2),'err2by2') % right y-axis
title('first select intensity threshold')
xlabel('spots');grid minor        

subplot(nFg,1,3); % PLOT #3
[hAx_spotSel1,hLineErr5,hLineIntPeak] = plotyy(s,err2G5by5_2plot(ixIn),s,intFit(ixIn)); hold on;
set(hLineErr5,'Color',[0 0 1],'MarkerEdgeColor','b','MarkerSize',mrkSz,'LineStyle','.','Marker','*')
set(hLineIntPeak,'Color',[0 0 0],'MarkerEdgeColor','r','MarkerSize',mrkSz,'LineStyle','.','Marker','*')        
legend('err2G5by5\_2','intFit');
ylabel(hAx_spotSel1(1),'err2G5by5\_2') % left y-axis
ylabel(hAx_spotSel1(2),'intFit') % right y-axis     
title('select a spot to display fitting')
xlabel('spots');grid minor        


subplot(nFg,1,4);  % PLOT #4
[hAx_spotSel2,hLineErr5_2,hLineIntPeak_2] = plotyy(s,err2G5by5_2plot(ixIn),s,intFit(ixIn)); hold on;
set(hLineErr5_2,'Color',[0 0 1],'MarkerEdgeColor','b','MarkerSize',mrkSz,'LineStyle','.','Marker','*')
set(hLineIntPeak_2,'Color',[0 0 0],'MarkerEdgeColor','r','MarkerSize',mrkSz,'LineStyle','.','Marker','*')        
legend('err2G5by5\_2','intFit');
ylabel(hAx_spotSel2(1),'err2G5by5\_2') % left y-axis
ylabel(hAx_spotSel2(2),'intFit') % right y-axis
title('select a spot to display fitting')
xlabel('spots');grid minor        

dispErr = 1;
subplot(nFg,1,5);  % PLOT #5
if dispErr
    [hLine1] = plot(s,isSpot(ixIn)); 
    %[hLine1] = plot(s,errCenter_(ixIn)); 
    
    set(hLine1,'Color',[0 0 1],'MarkerEdgeColor','k','MarkerSize',mrkSz,'LineStyle','.','Marker','*')
    legend('isSpot');
    xlabel('spots');grid minor
else
    [hAx,hLine1,hLine2] = plotyy(s,posCorr(ixIn),s,err2by2plot(ixIn)); hold on;
    set(hLine1,'Color',[0 0 1],'MarkerEdgeColor','k','MarkerSize',mrkSz,'LineStyle','.','Marker','*')
    set(hLine2,'Color',[0 0 0],'MarkerEdgeColor','r','MarkerSize',mrkSz,'LineStyle','.','Marker','*')        
    legend('posCorr','err2by2');
    ylabel(hAx(1),'posCorr') % left y-axis
    ylabel(hAx(2),'err2by2') % right y-axis
    xlabel('spots');grid minor
end
title('select a spot to display fitting')

subplot(nFg,1,6);  % PLOT #6
[hAx,hLine1,hLine2] = plotyy(s,err2G5by5_2plot(ixIn),s,err2by2plot(ixIn)); hold on;
set(hLine1,'Color',[0 0 1],'MarkerEdgeColor','b','MarkerSize',mrkSz,'LineStyle','.','Marker','*')
set(hLine2,'Color',[0 0 0],'MarkerEdgeColor','c','MarkerSize',mrkSz,'LineStyle','.','Marker','*')        
legend('err2G5by5\_2','err2by2');
ylabel(hAx(1),'err2G5by5\_2') % left y-axis
ylabel(hAx(2),'err2by2') % right y-axis
title('select a spot to display fitting')
xlabel('spots');grid minor

subplot(nFg,1,7);  % PLOT #7
[hAx,hLine1,hLine2] = plotyy(s,intFit(ixIn),s,intPeak(ixIn)); hold on;
set(hLine1,'Color',[0 0 1],'MarkerEdgeColor','r','MarkerSize',mrkSz,'LineStyle','.','Marker','*')
set(hLine2,'Color',[0 0 0],'MarkerEdgeColor','k','MarkerSize',mrkSz,'LineStyle','.','Marker','*')        
legend('intFit','intPeak');
ylabel(hAx(1),'intFit') % left y-axis
ylabel(hAx(2),'intPeak') % right y-axis
title('select a spot to display fitting')
xlabel('spots');grid minor


% select spot
[x_,y_]=ginput(1);
sp = round(x_);
isFit5by5Array = isFit5by5;

posFig = get(hAxErr(1),'Position');
x1Fig = posFig(1);
dxFig = posFig(3);

xSpot = x1Fig + sp*dxFig/numSpots;
hAnnotSpot = annotation('line',[xSpot xSpot],[0.1 0.95]);

if exist('fig1'), figure(fig1); else, fig1 = figure; end;
colormap('jet');
set(fig1,'Color',[0.2 0.2 0.2]*4)
axe=axes('Parent',fig1,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren','XLim',0.5+[0 m],'YLim',0.5+[0 n]);
ss = get(0,'MonitorPositions');set(fig1,'position',ss(1,:));        

q = 0; prnt = 0;
isErrThreshSelected = 0;
isSnap = 0;
while q == 0 
    xSpot = x1Fig + sp*dxFig/numSpots;
    set(hAnnotSpot,'X',[xSpot xSpot]);
    
    % mark the spot
    if sp < 0, sp =1; end;
    sp_ = sp2(sp);
    if isGenSpot, set(hText1,'BackgroundColor',[0 0 1]); end        % generated spot is detected

    if prnt nFigY=1; nFigX=4; close(fig1); fig1=figure; else nFigY=2;nFigX=4; end;
    %if ~isempty(isSpot), set(fig1,'Color',[0.2 0.2 0.2]);end; % nonfit
    tit = sprintf('5by5 Gaussian Fit');
    set(fig1,'Name',tit,'NumberTitle','off');       
    isFit5by5 = cell2mat(isFit5by5Array(sp));

    text1 = sprintf('spot:#%i,frame:%i,\n Ipeak:%.02f, X:%.02f,Y:%.02f'...
        ,sp,ixTrack(sp_),intPeak(sp),xBW(iyTrack(sp_),ixTrack(sp_)),yBW(iyTrack(sp_),ixTrack(sp_)));
    text2 = sprintf('PREFIT: err2G5by5:%.02f, err2by2:%.02f, \n ratio of intensities %.02f,dist:%.02f'...
        ,err2G5by5(sp),err2by2(sp),Irat(sp),dist(sp));
    text3 = sprintf('FIT: err2G5by5:%.02f, err2by2:%.02f, \n ratio of intensities %.02f,dist:%.02f'...
        ,err2G5by5_2(sp),err2by2_2(sp),Irat_2(sp),dist_2(sp));
    text4 = sprintf('isSpot:%i'...
        ,isSpot(sp));
    text5 = sprintf('exit flag:%s\n exit flag:%s\n exit flag:%s'...
        ,cell2mat(exitflag(sp,1)),cell2mat(exitflag(sp,2)),cell2mat(exitflag(sp,3)));
    text6 = sprintf('3by3 err:%s\n errCenter:%s '...
        ,isFit5by5,cell2mat(errCenter(sp)));
    %hText1 = uicontrol('style','text','BackgroundColor',[1 1 1],'String',text,'Position',[1300 500 360 80]);
    hText1 = uicontrol('style','text','BackgroundColor',[1 1 1],'String',text1,'Position',[1100 630 460 30]);
    hText2 = uicontrol('style','text','BackgroundColor',[1 1 1],'String',text2,'Position',[1100 590 460 30]);
    hText3 = uicontrol('style','text','BackgroundColor',[1 1 1],'String',text3,'Position',[1100 550 460 30]);
    hText4 = uicontrol('style','text','BackgroundColor',[1 1 1],'String',text4,'Position',[1100 520 460 20]);
    hText5 = uicontrol('style','text','BackgroundColor',[1 1 1],'String',text5,'Position',[1100 450 460 60]);
    hText6 = uicontrol('style','text','BackgroundColor',[1 1 1],'String',text6,'Position',[600 500 460 90]);
    figure(hFigCont); 
    hTextAct = uicontrol('style','text','BackgroundColor',[1 1 1],'String','select spot','Position',[100 1730 160 20]);
    hTextThr = uicontrol('style','text','BackgroundColor',[1 1 1],'String','threshold values:','Position',[400 90 160 20]);
    hTextInt = uicontrol('style','text','BackgroundColor',[1 1 1],'String','min int:','Position',[400 30 160 20]);
    hTextErr = uicontrol('style','text','BackgroundColor',[1 1 1],'String','max Err5:','Position',[400 60 160 20]);
    figure(fig1)
    
    %% FIG 1
    figIdatSHOW = subplot(nFigY,nFigX,1); 
    hIdatSHOW = imagesc(IdatSHOW(:,:,sp)); axis image;
    hold on; scatter(PX(sp)+2,PY(sp)+2); scatter(PX2by2(sp)+2,PY2by2(sp)+2,'x'); 
%    scatter(cc1disp(3,sp),cc1disp(5,sp),'+'); 
    hold off;
    colorbar;
    rectangle('Position', [2.5 2.5 5 5],'EdgeColor','r')
    if isempty(isFit5by5), rectangle('Position', [xCoor5by5(sp)-0.5 yCoor5by5(sp)-0.5 3 3],'EdgeColor','g','LineStyle','-.'); end;
%    if isempty(isFit5by5), rectangle('Position', [xCoor2by2(sp)-0.5 yCoor2by2(sp)-0.5 2 2],'EdgeColor','r','LineStyle','--'); end;
    title({'spot image';sprintf('fit2: err2G5by5:%.02f,err2by2:%.02f, x:%.02f ,y:%.02f',err2G5by5(sp),err2by2(sp),PX(sp),PY(sp))});

    %% SLIDER cross sections
    set(figIdatSHOW,'units','pixels');
    figPosIdatSHOW = get(figIdatSHOW,'Position');
    
    hCrossSectX = uicontrol('style','slider','units','pixel','Position',[200 figPosIdatSHOW(2)+35 20 figPosIdatSHOW(3)], 'SliderStep',[1/numSpots 10/numSpots]);
    hCrossSectY = uicontrol('style','slider','units','pixel','Position',[figPosIdatSHOW(1) 600 figPosIdatSHOW(3) 20], 'SliderStep',[1/numSpots 10/numSpots]);
    hXSectionButton = uicontrol('style','pushbutton','String','set sections','Position',[100 780 100 20],'Callback', 'updScrollState(hXSectionButton,hTextSpotNumState)');     
    
    ps = BigWindowSize; % section position
    posSect = (round([ps ps ps ps]/2));
    [sectDataX,sectDataY,sectFitX,sectFitY] = getCrossSect(IdatSHOW(:,:,sp),Ifit2SHOW_2(:,:,sp),posSect,1);
    subplot(nFigY,nFigX,2);  % 2 cross sections
    hSectDataX = plot(sectDataX);
    hold on; 
    hSectDataY = plot(sectDataY);
    hSectFitX = plot(sectFitX);
    hSectFitY = plot(sectFitY);
    hSect = [hSectDataX,hSectDataY,hSectFitX,hSectFitY];
    hold off;
    title({'cross sections';'__'});
    listenCrossSectionX = addlistener(hCrossSectX,'ActionEvent',@(hObject, event) updXsection(hObject, event,hCrossSectX,hCrossSectY,hSect,IdatSHOW(:,:,sp),Ifit2SHOW_2(:,:,sp),BigWindowSize,intMul));
    listenCrossSectionY = addlistener(hCrossSectY,'ActionEvent',@(hObject, event) updXsection(hObject, event,hCrossSectX,hCrossSectY,hSect,IdatSHOW(:,:,sp),Ifit2SHOW_2(:,:,sp),BigWindowSize,intMul));


    %% FIGS > 2
    subplot(nFigY,nFigX,3); imagesc(Ifit2SHOW_2(:,:,sp)); axis image; % 2 GAUS FIT
    hold on; scatter(cc2disp_2(8,sp),cc2disp_2(10,sp),'filled','b'); scatter(cc2disp_2(3,sp) ,cc2disp_2(5,sp),'filled','r'); hold off;
    colorbar;
    rectangle('Position', [2.5 2.5 5 5],'EdgeColor','r')
    if isempty(isFit5by5), rectangle('Position', [xCoor5by5(sp)-0.5 yCoor5by5(sp)-0.5 3 3],'EdgeColor','g','LineStyle','-.'); end;
%    if isempty(isFit5by5), rectangle('Position', [xCoor2by2(sp)-0.5 yCoor2by2(sp)-0.5 2 2],'EdgeColor','r','LineStyle','--'); end;
    title({'2GausFit w\ constraints'; sprintf('gaus1(red): I:%.02f, x:%.02f ,y:%.02f, sgx:%.02f,sgy:%.02f',cc2disp_2(2,sp),cc2disp_2(3,sp),cc2disp_2(5,sp),cc2disp_2(4,sp),cc2disp_2(6,sp));sprintf('gaus2(blue):  I:%.02f, x:%.02f ,y:%.02f, sgx:%.02f,sgy:%.02f',cc2disp_2(7,sp),cc2disp_2(8,sp),cc2disp_2(10,sp),cc2disp_2(9,sp),cc2disp_2(11,sp))});
    if ~prnt % display mode
%         subplot(nFigY,nFigX,4); imagesc(Ifit1SHOW(:,:,sp)); axis image; % 1 GAUS FIT
%         hold on; scatter(cc1disp(3,sp),cc1disp(5,sp),'filled'); hold off;
%         colorbar;
%         rectangle('Position', [2.5 2.5 5 5],'EdgeColor','r')
%         if isempty(isFit5by5), rectangle('Position', [xCoor5by5(sp)-0.5 yCoor5by5(sp)-0.5 3 3],'EdgeColor','g','LineStyle','-.'); end;
%         if isempty(isFit5by5), rectangle('Position', [xCoor2by2(sp)-0.5 yCoor2by2(sp)-0.5 2 2],'EdgeColor','r','LineStyle','--'); end;
%         title(sprintf('gaus:  I:%.02f, x:%.02f ,y:%.02f, sgx:%.02f,sgy:%.02f',cc1disp(2,sp),cc1disp(3,sp),cc1disp(5,sp),cc1disp(4,sp),cc1disp(6,sp)));

        subplot(nFigY,nFigX,5); 
        imgDiff = (IdatSHOW(:,:,sp) - Ifit2SHOW_2(:,:,sp));
        %mnImg = min(imgDiff(:)); mxImg = max(imgDiff(:));
        %imgDiff = (imgDiff - mnImg)/(mxImg-mnImg)*63+1;
        imagesc(imgDiff(y5y1:y5y2,x5x1:x5x2));
        axis image
        hC = colorbar; 
        %[ticks64, tickPos, ticks ] = plotColorBar(imgDiff,6,0);
        %set(hC,'Ytick',tickPos,'YTicklabel',ticks);

        % frame
        subplot(nFigY,nFigX,6); 
        binnedFrm = ixTrack(sp_)-1;
        for iAv = 1 : binFrame
            frmRead = binnedFrm*binFrame+iAv;
            temp(:,:,iAv) = double(imread(fname,frmRead)); 
        end        
        img1 = mean(temp(yy1:yy2,xx1:xx2,:),3);
        hImg = imagesc(img1); % axis image; 
        axis image; hold on; scatter(xBW(iyTrack(sp_),ixTrack(sp_)),yBW(iyTrack(sp_),ixTrack(sp_)),'w'); hold off

        % neighbouring frames
        subplot(nFigY,nFigX,[7 8]); 
        if ixTrack(sp_) > 1
            for iAv = 1 : binFrame
                frmRead = (binnedFrm-1)*binFrame+iAv;
                temp(:,:,iAv) = double(imread(fname,frmRead)); 
            end        
            img2 = mean(temp(yy1:yy2,xx1:xx2,:),3);        
        else
            img2=[]; 
        end
        
        for iAv = 1 : binFrame
            frmRead = (binnedFrm+1)*binFrame+iAv;
            temp(:,:,iAv) = double(imread(fname,frmRead)); 
        end        
        img3 = mean(temp(yy1:yy2,xx1:xx2,:),3);        
        
        IMG = [img2 img1 img3];
        hImg = imagesc(IMG); %axis image; 
        [szY,szX]=size(IMG);
        
        if ixTrack(sp_) > 1
            hold on; 
            scatter([szX/6 szX/2 szX*5/6]+0.5,[szY szY szY]/2+0.5,'r*');
            hold off;
        end
        axis image; 
        posLine1 = size(img1,2)+0.5;
        lenLine =  size(img1,1);
        posLine2 = size(img1,2)*2+0.5;
        hLine1=line([posLine1 posLine1],[0.5 0.5+lenLine],'Color',[0 0 0]);
        hLine2=line([posLine2 posLine2],[0.5 0.5+lenLine],'Color',[0 0 0]);
        
    else
        % update figure info
        set(fig1,'position',[1 1 1800 600]); 
        set(hText1,'Position',get(hText1,'Position')+[200 -200 0 0])
        set(hText2,'Position',get(hText2,'Position')+[200 -200 0 0])
        set(hText3,'Position',get(hText3,'Position')+[200 -200 0 0])
        set(hText4,'Position',get(hText4,'Position')+[200 -200 0 0])
        set(hText5,'Position',get(hText5,'Position')+[200 -200 0 0])
        set(hText6,'Position',get(hText6,'Position')+[200 -200 0 0])
        % print images
        imgFig = getframe(gcf,[200,150,1600,350]);
        if sp == 1
            imwrite(imgFig.cdata,'spotFitView.tif') 
        else
            imwrite(imgFig.cdata,'spotFitView.tif','WriteMode','append') 
        end
    end
    replot = 0;
    if prnt, replot = 1; sp = sp + 1; end;
    if sp==SP, q = 1; end % combine frames in a stack
    
    while replot == 0
        btn = 0;
        replot = 1;
        figure(hFigCont); set(hTextAct,'String','press a to select or q to quit')
        textState = get(hTextSpotNumState,'String'); % 
        if strcmp(textState,'select spot')
            while btn == 0
                btn = waitforbuttonpress;
                k = get(hFigCont,'CurrentCharacter');
            end
            switch lower(k)
                case 'a' % select another spot
                    figure(hFigCont)
                    [x_,y_]=ginput(1);
                    sp = round(x_);
                case 'q' % quit
                    set(hTextAct,'String','EXITING ...')
                    q=1;
                case 'p' % print spots
                    set(hTextAct,'String','PRINTING ...')
                    sp = 1;
                    prnt=1;
            end            
        elseif strcmp(textState,'scroll')
            while sp == spPrev
                slideVal = get(hSpotSlider,'Value');
                sp = round((numSpots-1)*slideVal)+1;
                pause(0.1)
            end
            slideVal = get(hSpotSlider,'Value');
            sp = round((numSpots-1)*slideVal)+1;
            spPrev = sp;
            figure(fig1)
            break;            
        elseif strcmp(textState,'set cross section')
            figure(fig1)
            break
            
        end
        if (gca == hAx_spotSel1(1) || gca == hAx_spotSel1(2)) 
            
            if isErrThreshSelected
                ixDisp = ixSel;
            else
                ixDisp = ixInt5_2;
            end
            isSnap = 1;     
        elseif (gca == hAx_spotSel2(1) || gca == hAx_spotSel2(2)) 
            if isErrThreshSelected
                ixDisp = ixNonSel; % discarded
            else
                ixDisp = ixInt5_2;
            end
            isSnap = 1;
        else
            isSnap = 0;
        end
        if isSnap
            [valMin ixMin] = min(abs(ixDisp - x_ ));
            sp = ixDisp(ixMin);
        end


        while y_ < 1 && y_ > 0 && (gca == hAxInt(1) || gca == hAxInt(2))     % intensity                 
            title(hAx_spotSel1(1),'spots below the intensity threshold, select error threshold to further filter')
            title(hAx_spotSel2(1),'spots below the intensity threshold, select error threshold to further filter')            
            set(hTextAct,'String','select intensity threshold')
            err5 = err2G5by5_2plot;
            ixInt5 = find(intFit>=y_); % selected
            ixInt5_2 = find(intFit<y_); % some will be selected by Err performance
            set(hLineErr5,      'YData',err5(ixInt5_2));
            set(hLineErr5,      'XData',ixInt5_2);
            set(hLineIntPeak,   'YData',intFit(ixInt5_2));
            set(hLineIntPeak,   'XData',ixInt5_2);
            set(hLineErr5_2,      'YData',err5(ixInt5_2));
            set(hLineErr5_2,      'XData',ixInt5_2);
            set(hLineIntPeak_2,   'YData',intFit(ixInt5_2));
            set(hLineIntPeak_2,   'XData',ixInt5_2);
            if exist('hline2')
                set(hline2,'YData',[y_ y_])
            else
                hline2 = line([0,s(end)],[y_ y_]);
            end
            intThreshold = y_;% *mxIntPeak
            set(hTextInt,'String',sprintf('min int:%.02f',y_));                    
            [x_,y_]=ginput(1);
            replot = 0;
            %title(hAxErr(1),'spots below intensity threshold, select error threshold to further filter')

        end
        while y_ < 1 && y_ > 0 && (gca == hAxErr(1) || gca == hAxErr(2))
            title(hAx_spotSel1(1),'non selected spots')
            title(hAx_spotSel2(1),'spots selected to be discarded as non-fit')
            set(hTextAct,'String','select error threshold')
            err5 = err2G5by5_2plot(ixInt5_2);
            ixErr5 = find(err5<=y_); % also selected
            ixErr5_2 = find(err5>y_); % not selected
            ixSel = sort([ixInt5_2(ixErr5) ixInt5]);
            ixNonSel = ixInt5_2(ixErr5_2);
            %ixNonSel = sort([ixErr5_2 ixInt5_2]);

            set(hLineErr5,      'YData',err2G5by5_2plot(ixSel));
            set(hLineErr5,      'XData',ixSel);
            set(hLineIntPeak,   'YData',intFit(ixSel));
            set(hLineIntPeak,   'XData',ixSel);
            set(hLineErr5_2,      'YData',err2G5by5_2plot(ixNonSel));
            set(hLineErr5_2,      'XData',ixNonSel);
            set(hLineIntPeak_2,   'YData',intFit(ixNonSel));
            set(hLineIntPeak_2,   'XData',ixNonSel);
            if exist('hline')
                set(hline,'YData',[y_ y_])
            else
                hline = line([0,s(end)],[y_ y_]);
            end
            set(hTextAct,'String','select error/intensity threshold')
            errThreshold = y_*(mxErr-mnErr)+mnErr;
            set(hTextErr,'String',sprintf('max Err5:%.02f',y_));
            [x_,y_]=ginput(1);
            % set(hFigCont,'Position', figPos);        
            replot = 0;
            isErrThreshSelected = 1;
        end
        if exist('ixSel')
            subplot(nFg,1,5);
            hold on; 
            pp = ones(size(ixSel))+0.2;
            plot(ixSel,pp,'r-v','MarkerSize',3);
            legend('isSpot','passed');
            ylim([0.2 2]);
            hold off;   
        end
        if q ~= 1, set(hTextAct,'String','select spot'); end;
        figure(fig1)
    end 
end
save('spotSelVal','ixSel','ixNonSel','errThreshold','intThreshold');
clear;
load spotSelVal; load spotSel; save spotSel; 
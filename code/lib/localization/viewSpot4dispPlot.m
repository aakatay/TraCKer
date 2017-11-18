
% data
load spotInfo;
numFrmDisp = 1000;
numFrmDisp = 200;
cropDisp

gry = 0.6; % text background color coeff

y7y1=x7x1;y7y2=x7x2;
m = 7; n = 7; m2 = m + 2; n2 = n + 2;

px0 = double(ceil(m/2));
% define output XY array
[iyTrack, ixTrack] = find(xBW>0); % positions of spot data in the X array

frstSpot = 1;
ixIn = frstSpot:SP;
spPrev = 0;

% inactive options
isGenSpot = 0;

numSpots = size(fitVal,2);

for i = 1 : length(errCenter)
    if isempty(cell2mat(errCenter(i)))
    errCenter_(i) = nan;
    else
    errCenter_(i) = 1;        
    end
end

%intFitNorm(200)=0;
%err2by2plot(200)=0;

tit2 = 'controls';
hFigCont = figure('DoubleBuffer','on','Menubar','none','Name',tit2,'NumberTitle','off','Colormap',gray(256));
figPos = [1921 -461 1080 1844];
set(hFigCont,'Position', figPos);        

% error plots
nFg = 7;
errFit7by7(errFit7by7==0) = nan;
hText_ = uicontrol('style','text','BackgroundColor',[1 1 1]*gry,'String',sprintf('max. intPeak:%.02f',max(intPeak(:))),'Position',[900 1730 160 20]);

% normalize
intPeakMax = max(intPeak(:));
intPeakNorm = intPeak/intPeakMax;

intFit = fitVal(2,:);
intFitMax = max(intFit);
intFitNorm = intFit / intFitMax;


%mn = min(err2by2(:)); mx = max(err2by2(:)); err2by2plot = (err2by2-mn)/(mx-mn);
err2by2plot = err2by2;
%hText_ = uicontrol('style','text','BackgroundColor',[1 1 1]*gry,'String',sprintf('err2by2: %.02f - %.02f',mn,mx),'Position',[100 1500 160 20]);        

mn = min(errFit7by7(:)); mx = max(errFit7by7(:)); errFit7by7plot = (errFit7by7-mn)/(mx-mn);
mnErr = mn; mxErr = mx; 

hText_ = uicontrol('style','text','BackgroundColor',[1 1 1]*gry,'String',sprintf('errFit7by7: %.02f - %.02f',mn,mx),'Position',[800 1500 160 20]);        

hScrollButton = uicontrol('style','pushbutton','String','start scroll','Position',[100 1780 100 20],'Callback', 'updScrollState(hScrollButton,hTextSpotNumState)');     
hTextSpotNumState = uicontrol('style','text','BackgroundColor',[1 1 1]*gry,'String',sprintf('select spot'),'Position',[250 1780 100 20]);     
hTextSpotNum = uicontrol('style','text','BackgroundColor',[1 1 1]*gry,'String',sprintf('spot#: %i',0),'Position',[400 1780 100 20]);     
hSpotSlider = uicontrol('style','slider','units','pixel','Position',[100 1750 900 20], 'SliderStep',[1/numSpots 10/numSpots]);
listenSpot = addlistener(hSpotSlider,'ActionEvent',@(hObject, event) updSpotNum(hObject, event,hTextSpotNum,numSpots));
%listenScrollBtn = addlistener(hScrollButton,'ActionEvent',@(hObject, event) updScrollState(hObject, event,hScrollButton,hTextSpotNumState));


%% spot selection figures
mrkSz = 2;
subplot(nFg,1,1);   s = ixIn; % PLOT #1
[hAxInt,hLine1,hLine2] = plotyy(s,intPeakNorm(ixIn),s,intFitNorm(ixIn)); hold on;
set(hLine1,'Color','w','MarkerSize',mrkSz,'LineStyle','.','Marker','*')
set(hLine2,'Color','r','MarkerSize',mrkSz,'LineStyle','.','Marker','*')
legend('intPeakNorm','intFitNorm');
ylim([0 1])
ylabel(hAxInt(1),'intPeakNorm') % left y-axis
ylabel(hAxInt(2),'intFitNorm') % right y-axis
title('select intensity threshold to filter spots')
xlabel('spots');grid minor
setFigColor(gcf,hAxInt)

subplot(nFg,1,2); % PLOT #2
[hAxErr,hLine1,hLine2] = plotyy(s,errFit7by7plot(ixIn),s,err2by2plot(ixIn)); hold on;
%plot(s,errFit7by7);
set(hLine1,'Color','c','MarkerSize',mrkSz,'LineStyle','.','Marker','*')
set(hLine2,'Color','g','MarkerFaceColor',[0.5 0.5 1],'MarkerSize',mrkSz,'LineStyle','.','Marker','*')        
legend('errFit7by7','err2by2');
ylabel(hAxErr(1),'errFit7by7') % left y-axis
ylabel(hAxErr(2),'err2by2') % right y-axis
title('first select intensity threshold')
xlabel('spots');grid minor        
setFigColor(gcf,hAxErr)

subplot(nFg,1,3); % PLOT #3
[hAx_spotSel1,hLineErr5,hLineIntPeak] = plotyy(s,errFit7by7plot(ixIn),s,intFitNorm(ixIn)); hold on;
set(hLineErr5,'Color','c','MarkerSize',mrkSz,'LineStyle','.','Marker','*')
set(hLineIntPeak,'Color','r','MarkerSize',mrkSz,'LineStyle','.','Marker','*')        
legend('errFit7by7','intFitNorm');
ylabel(hAx_spotSel1(1),'errFit7by7') % left y-axis
ylabel(hAx_spotSel1(2),'intFitNorm') % right y-axis     
title('select a spot to display fitting')
xlabel('spots');grid minor        
setFigColor(gcf,hAx_spotSel1)


subplot(nFg,1,4);  % PLOT #4
[hAx_spotSel2,hLineErr5_2,hLineIntPeak_2] = plotyy(s,errFit7by7plot(ixIn),s,intFitNorm(ixIn)); hold on;
set(hLineErr5_2,'Color','c','MarkerSize',mrkSz,'LineStyle','.','Marker','*')
set(hLineIntPeak_2,'Color','r','MarkerSize',mrkSz,'LineStyle','.','Marker','*')        
legend('errFit7by7','intFitNorm');
ylabel(hAx_spotSel2(1),'errFit7by7') % left y-axis
ylabel(hAx_spotSel2(2),'intFitNorm') % right y-axis
title('select a spot to display fitting')
xlabel('spots');grid minor        
setFigColor(gcf,hAx_spotSel2);

dispErr = 1;
subplot(nFg,1,5);  % PLOT #5
if dispErr
    [hLine1] = plot(s,ones(size(ixIn))); 
    %[hLine1] = plot(s,errCenter_(ixIn)); 
    
    set(hLine1,'Color','w','MarkerSize',mrkSz,'LineStyle','.','Marker','*')
    legend('isDoubleSpot');
    xlabel('spots');
    ylim([0 2])
    setFigColor(gcf,gca)
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
[hAx,hLine1,hLine2] = plotyy(s,errFit7by7plot(ixIn),s,err2by2plot(ixIn)); hold on;
set(hLine1,'Color','c','MarkerSize',mrkSz,'LineStyle','.','Marker','*')
set(hLine2,'Color','g','MarkerSize',mrkSz,'LineStyle','.','Marker','*')        
legend('errFit7by7','err2by2');
ylabel(hAx(1),'errFit7by7') % left y-axis
ylabel(hAx(2),'err2by2') % right y-axis
title('select a spot to display fitting')
xlabel('spots');grid minor
setFigColor(gcf,hAx)

subplot(nFg,1,7);  % PLOT #7
[hAx,hLine1,hLine2] = plotyy(s,fitVal(1,ixIn),s,IntFlatBckgrnd(ixIn)); hold on;
set(hLine1,'Color',[0.4 0.4 0],'MarkerSize',mrkSz,'LineStyle','.','Marker','*')
set(hLine2,'Color',[0 0.4 0.4 ],'MarkerSize',mrkSz,'LineStyle','.','Marker','*')        
legend('bckgrndFit','bckgrndMedian');
ylabel(hAx(1),'bckgrndFit') % left y-axis
ylabel(hAx(2),'bckgrndMedian') % right y-axis
title('select a spot to display fitting')
xlabel('spots');grid minor
setFigColor(gcf,hAx)


%% select spot
[x_,y_]=ginput(1);
sp = round(x_);
 
posFig = get(hAxErr(1),'Position');
x1Fig = posFig(1);
dxFig = posFig(3);

xSpot = x1Fig + sp*dxFig/numSpots;
hAnnotSpot = annotation('line',[xSpot xSpot],[0.1 0.95]);
set(hAnnotSpot,'Color',[0 0 0.6])

if exist('fig1'), figure(fig1); else fig1 = figure; end;
colormap('jet');
set(fig1,'Color',[0.2 0.2 0.2]*4)
axe=axes('Parent',fig1,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren','XLim',0.5+[0 m],'YLim',0.5+[0 n]);
ss = get(0,'MonitorPositions');set(fig1,'position',ss(1,:));        

q = 0; prnt = 0;
isErrThreshSelected = 0;
isSnap = 0;
countSp = 1;
while q == 0 
    xSpot = x1Fig + sp*dxFig/numSpots;
    set(hAnnotSpot,'X',[xSpot xSpot]);
    
    
    if sp < 0, sp =1; end;
    % mark the spot
    if 0 
        switch rem(countSp,3)
            case 1
                sp = 47;
            case 2
                sp = 41;
            case 0
                sp = 42;
        end
        countSp = countSp + 1;
    end
%79%60; % ===============================================================
%27
%79
    if isGenSpot, set(hText1,'BackgroundColor',[0 0 1]); end        % generated spot is detected

    %% displayed images
    fitValExt           = [0 0 1 0 1 0 0 1 0 1 0 0 1 0 1 0 0 1 0 1 0 0 1 0 1 0 0 1 0 1 0 0 1 0 1 0 ]'; % coor extension
    nSptsReg = nSpts(1,sp);
    nSptsNonReg = nSpts(2,sp);
    regDispOffset       = [0 0 1 0 1 0 0 1 0 1 0 0 1 0 1 0 0 1 0 1 0 0 1 0 1 0 0 1 0 1 0]';
    nonRegDispOffset    = [0 0 1 0 1 0 0 1 0 1 0 0 1 0 1 0 0 1 0 1 0 0 1 0 1 0 0 1 0 1 0 0 1 0 1 0 ]';
    if isempty(nSptsNonReg)
        numParam = 1+5;
    else
        numParam = 1+5+5*nSptsNonReg;
    end
    xPosVec = [3,8,13,18,23,28,33,38];
    xPosVec = xPosVec(xPosVec<=numParam);
    
    nonRegDispOffset = nonRegDispOffset(1:numParam);
    fitValDisp =  fitVal(1:numParam,sp)+nonRegDispOffset;
    fitValDisp(fitVal(1:numParam,sp)==0) = 0;
    ftX = fitValDisp(xPosVec);    % fit positions
    ftY = fitValDisp(xPosVec+2);
    G0regDisp = G0reg(:,sp)+regDispOffset;
    G0nonRegDisp = G0nonReg(1:numParam,sp)+nonRegDispOffset;
    ftX0 = G0nonRegDisp(xPosVec); % fit starting positions
    ftY0 = G0nonRegDisp(xPosVec+2); 
    
    % fit function: funSpNonReg
    load('spotSel.mat', 'intP')
    if ~exist('txFitSpNonRegFunc')
        [funText, txFitSpNonRegFunc] = genFitFunc(nSptsNonReg); % function text
    else
        [funText, ~] = genFitFunc(nSptsNonReg,txFitSpNonRegFunc); % function text
    end
    eval(funText);

    regX = G0regDisp([3,8,13,18,23,28]); % pre. registered spots
    regY = G0regDisp([3,8,13,18,23,28]+2);
    regX = regX(regX>1);
    regY = regY(regY>1);
    
    % 1- fit : single fluo.
    fitPeak = fitValDisp;
    fitPeak([1,7:end])= 0;
    Ifit2SHOWpeak_ = funSpNonReg(fitPeak,xSHOW);
    Ifit2SHOWpeak=reshape(Ifit2SHOWpeak_,[n2 m2]);%gaussian reshaped as matrix
    Ifit2SHOWpeak = Ifit2SHOWpeak*intPeak(sp);

    % 2- fit : single fluo. with background
    Ifit2SHOW_=funSpNonReg(fitValDisp,xSHOW); %your fitted gaussian in vector
    Ifit2SHOW=reshape(Ifit2SHOW_,[n2 m2]); %gaussian reshaped as matrix
    Ifit2SHOW = (Ifit2SHOW+fit9by9reg(:,:,sp)+IntFlatBckgrndDiff9by9(:,:,sp)...
        +IntFlatBckgrnd(sp)+data9by9blur(:,:,sp))*intPeak(sp);

    if prnt nFigY=1; nFigX=4; close(fig1); fig1=figure; else nFigY=2;nFigX=4; end;
    %if ~isempty(isSpot), set(fig1,'Color',[0.2 0.2 0.2]);end; % nonfit
    tit = sprintf('5by5 Gaussian Fit');
    set(fig1,'Name',tit,'NumberTitle','off');    


        
    %% texts in the figure
    text1 = sprintf('spot:#%i,frame:%i,\n Ipeak:%.02f, X:%.02f,Y:%.02f'...
        ,sp,frm(sp),intPeakNorm(sp),xBW(sp),yBW(sp));
    text2 = sprintf('FIT: errFit7by7:%.02f, err2by2:%.02f'...
        ,errFit7by7(sp),err2by2(sp));
    text3Line1 = sprintf('bounds:  __sp1__     __sp2__     __sp3__ ');
    text3Line2 = sprintf('posX:  %.02f-%.02f,  %.02f-%.02f,  %.02f-%.02f',LB(3,sp),UB(3,sp),LB(8,sp),UB(8,sp),LB(13,sp),UB(13,sp));
    text3Line3 = sprintf('posY:  %.02f-%.02f,  %.02f-%.02f,  %.02f-%.02f',LB(5,sp),UB(5,sp),LB(10,sp),UB(10,sp),LB(15,sp),UB(15,sp));
    text3Line4 = sprintf('sigX:  %.02f-%.02f,  %.02f-%.02f,  %.02f-%.02f',LB(4,sp),UB(4,sp),LB(9,sp),UB(9,sp),LB(14,sp),UB(14,sp));
    text3Line4 = sprintf('sigY:  %.02f-%.02f,  %.02f-%.02f,  %.02f-%.02f',LB(6,sp),UB(6,sp),LB(11,sp),UB(11,sp),LB(16,sp),UB(16,sp));
    text3Line5 = sprintf('INT:   %.02f-%.02f,  %.02f-%.02f,  %.02f-%.02f',LB(2,sp),UB(2,sp),LB(7,sp),UB(7,sp),LB(12,sp),UB(12,sp));
    
    text3 = sprintf('%s \n %s \n %s \n %s \n %s \n',text3Line1,text3Line2,text3Line3,text3Line4,text3Line5);
    text4 = sprintf('nSptsNonReg: %i,nSptsReg: %i,',nSptsNonReg,nSptsReg);
    text5 = sprintf('exit flag:%s\n ',cell2mat(exitflag(sp)));
    text6 = sprintf('3by3 err:__\n errCenter:%s '...
        ,cell2mat(errCenter(sp)));
    del = 50;
    %hText1 = uicontrol('style','text','BackgroundColor',[1 1 1]*gry,'String',text,'Position',[1300 500 360 80]);
    xpos = 1500; 
    wh = 250;
    hText1 = uicontrol('style','text','BackgroundColor',[1 1 1]*gry,'String',text1,'Position',[xpos 630 wh 30]);
    hText2 = uicontrol('style','text','BackgroundColor',[1 1 1]*gry,'String',text2,'Position',[xpos 590 wh 30]);
    hText3 = uicontrol('style','text','BackgroundColor',[1 1 1]*gry,'String',text3,'Position',[xpos 550-del wh 30+del]);
    hText4 = uicontrol('style','text','BackgroundColor',[1 1 1]*gry,'String',text4,'Position',[xpos 520-del wh 20]);
    hText5 = uicontrol('style','text','BackgroundColor',[1 1 1]*gry,'String',text5,'Position',[xpos 450-del wh 20]);
    hText6 = uicontrol('style','text','BackgroundColor',[1 1 1]*gry,'String',text6,'Position',[1200 500 260 60]);
    
    % TEXT: spot positions
    LB_ = LB(:,sp)+regDispOffset(1:size(LB,1));
    UB_ = UB(:,sp)+regDispOffset(1:size(LB,1));
    nS = 5;
    viewSpot4dispPlot_text1; % CALLS : displays text for the fitted spots margins and fit values
    
    %% draw arrows between the starting and fit position
    
    
    % reg. spots
    clear x1spt y1spt x2spt y2spt
    for i = 1:nSptsReg % pre-fit spots
        x1spt(i) = double(xBW(spReg7by7(i,sp))+px0-xBW(sp)+1);
        y1spt(i) = double(yBW(spReg7by7(i,sp))+px0-yBW(sp)+1);
    end
    x2spt = regX;
    y2spt = regY;

    % fit spots
    x2spt(nSptsReg+1:nSptsReg+nSptsNonReg+1) = ftX;
    y2spt(nSptsReg+1:nSptsReg+nSptsNonReg+1) = ftY;
    x1spt(nSptsReg+1:nSptsReg+nSptsNonReg+1) = ftX0;
    y1spt(nSptsReg+1:nSptsReg+nSptsNonReg+1) = ftY0;
    
%% DEBUG:
if 0
    % main spot
    [xBW(sp) yBW(sp)]
    [Xgaus(sp) Ygaus(sp)]

    % 1st neigh spot
    spReg7by7(1,sp)
    [xBW(spReg7by7(1,sp)) yBW(spReg7by7(1,sp))]
    [Xgaus(spReg7by7(1,sp)) Ygaus(spReg7by7(1,sp))]
    G0reg([3,5],sp)+double([xBW(spReg7by7(1,sp)) yBW(spReg7by7(1,sp))]'-px0) 
    
    [fitVal(8,spReg) fitVal(10,spReg)]
    
   
    [x1spt' x2spt' y1spt' y2spt']
    
    G0reg(:,sp)
end    
    

    figure(hFigCont); 
    hTextAct = uicontrol('style','text','BackgroundColor',[1 1 1]*gry,'String','select spot','Position',[100 1730 160 20]);
    hTextThr = uicontrol('style','text','BackgroundColor',[1 1 1]*gry,'String','threshold values:','Position',[400 90 160 20]);
    hTextInt = uicontrol('style','text','BackgroundColor',[1 1 1]*gry,'String','min int:','Position',[400 30 160 20]);
    hTextErr = uicontrol('style','text','BackgroundColor',[1 1 1]*gry,'String','max Err5:','Position',[400 60 160 20]);
    
    figure(fig1)
    
    %% spot figures
    %% FIG 1
    figIdatSHOW = subplot(nFigY,nFigX,1); 
    hIdatSHOW = imagesc(spotWin(:,:,sp)); axis image;
    gcaImg(1) = gca;
    set(gca,'XColor',[1 1 1],'YColor',[1 1 1])
    hold on; 
    for i = 1: size(x1spt,2) % draw lines betw. starting and fit values
        line([x1spt(i) x2spt(i)],[y1spt(i) y2spt(i)]);
    end
    
    % 
    dxTx = 0.4; % text position offset
    txSptsNonReg = ['m';'1';'2';'3';'4';'5';'6'];
    scatter(ftX,ftY,'k.');
    htnonreg = text(ftX+dxTx,ftY+dxTx,txSptsNonReg(1:nSptsNonReg+1));
    set(htnonreg,'Color',[0 0 0],'FontWeight','bold','HorizontalAlignment','center')
    
    dxTx = -dxTx; 
    txSptsReg = ['1';'2';'3';'4';'5';'6'];    
    scatter(regX,regY,'k.');
    htreg = text(regX+dxTx,regY+dxTx,txSptsReg(1:nSptsReg),'EdgeColor',[0 0 0]);
    set(htreg,'Color',[0.5 0 0.5],'FontWeight','bold','HorizontalAlignment','center');
    
    %scatter(PX2(sp)+2,PY2(sp)+2,'x'); scatter(PX3(sp)+2,PY3(sp)+2,'cx'); scatter(cc1disp(3),cc1disp(5),'+'); 
    hold off;

    colorbar;
    rectangle('Position', [2.5 2.5 5 5],'EdgeColor','r')
    rectangle('Position', [xCoor3by3(sp)-0.5 yCoor3by3(sp)-0.5 3 3],'EdgeColor','g','LineStyle','-.'); 
    rectangle('Position', [xCoor2by2(sp)-0.5 yCoor2by2(sp)-0.5 2 2],'EdgeColor','r','LineStyle','--'); 
    %title({'spot image';sprintf('fit1: err1G5by5:%.02f, x:%.02f ,y:%.02f',err1G5by5(sp),cc1(3,sp),cc1(5,sp));sprintf('fit2: errFit7by7:%.02f,err2by2:%.02f, x:%.02f ,y:%.02f',errFit7by7(sp),err2by2(sp),PX(sp),PY(sp))});
    htt = title('image');
    set(htt,'Color',[1 1 1]);
    % SLIDER cross sections
    set(figIdatSHOW,'units','pixels');
    figPosIdatSHOW = get(figIdatSHOW,'Position');
    
    hCrossSectX = uicontrol('style','slider','units','pixel','Position',[200 figPosIdatSHOW(2)+65 20 figPosIdatSHOW(3)], 'SliderStep',[1/numSpots 10/numSpots]);
    hCrossSectY = uicontrol('style','slider','units','pixel','Position',[figPosIdatSHOW(1) 600+65 figPosIdatSHOW(3) 20], 'SliderStep',[1/numSpots 10/numSpots]);
    hXSectionButton = uicontrol('style','pushbutton','String','set sections','Position',[100 780 100 20],'Callback', 'updScrollState(hXSectionButton,hTextSpotNumState)');     
    
    %% FIG 2
    subplot(nFigY,nFigX,2);  % 2 cross sections
    ps = BigWindowSize; % section position
    posSect = (round([ps ps ps ps]/2));
   
    [sectDataX,sectDataY,sectFitX,sectFitY] = getCrossSect(spotWin(:,:,sp),Ifit2SHOW,posSect,1);
    hSectDataX = plot(sectDataX,'r');
    hAx = gca;
    setFigColor(gcf,hAx);
    hold on; 
    hSectDataY = plot(sectDataY,'Color',[1 0.5 0.5]);
    hSectFitX = plot(sectFitX);
    hSectFitY = plot(sectFitY,'c');
    hSect = [hSectDataX,hSectDataY,hSectFitX,hSectFitY];
    hold off;
    htt = title({'cross sections';'__'});
    set(htt,'Color',[1 1 1]);
    hLegend = legend('dataX','dataY','fitX','fitY');
    set(hLegend,'TextColor',[1 1 1]);
    %%
    listenCrossSectionX = addlistener(hCrossSectX,'ActionEvent',@(hObject, event) updXsection(hObject, event,hCrossSectX,hCrossSectY,hSect,spotWin(:,:,sp),Ifit2SHOW,BigWindowSize,intPeak));
    listenCrossSectionY = addlistener(hCrossSectY,'ActionEvent',@(hObject, event) updXsection(hObject, event,hCrossSectX,hCrossSectY,hSect,spotWin(:,:,sp),Ifit2SHOW,BigWindowSize,intPeak));


    %% FIGS > 2
    subplot(nFigY,nFigX,3); imagesc(Ifit2SHOW); axis image; % 2 GAUS FIT
    set(gca,'XColor',[1 1 1],'YColor',[1 1 1])
    gcaImg(2) = gca;
    hold on; 
    %scatter(fitValDisp(3),fitValDisp(10),'filled','b'); 
    scatter(fitValDisp(3) ,fitValDisp(5),'filled','r'); hold off;
    colorbar;
    rectangle('Position', [2.5 2.5 5 5],'EdgeColor','r')
    rectangle('Position', [xCoor3by3(sp)-0.5 yCoor3by3(sp)-0.5 3 3],'EdgeColor','g','LineStyle','-.'); 
    rectangle('Position', [xCoor2by2(sp)-0.5 yCoor2by2(sp)-0.5 2 2],'EdgeColor','r','LineStyle','--'); 
    htt = title({'2GausFit w\ constraints'; sprintf('gaus1(red): I:%.02f, x:%.02f ,y:%.02f, sgx:%.02f,sgy:%.02f',fitValDisp(2),fitValDisp(3),fitValDisp(5),fitValDisp(4),fitValDisp(6));...
        %sprintf('gaus2(blue):  I:%.02f, x:%.02f ,y:%.02f, sgx:%.02f,sgy:%.02f',fitValDisp(7),fitValDisp(8),fitValDisp(10),fitValDisp(9),fitValDisp(11))
        });
    set(htt,'Color',[1 1 1]);
    %% FIG 4
    if ~prnt % display mode
        subplot(nFigY,nFigX,4); imagesc(Ifit2SHOWpeak); axis image; % 1 GAUS FIT
        set(gca,'XColor',[1 1 1],'YColor',[1 1 1])
        gcaImg(3) = gca;
        %hold on; scatter(cc1disp(3),cc1disp(5),'filled'); hold off;
        colorbar;
        rectangle('Position', [2.5 2.5 5 5],'EdgeColor','r')
        rectangle('Position', [xCoor3by3(sp)-0.5 yCoor3by3(sp)-0.5 3 3],'EdgeColor','g','LineStyle','-.'); 
        rectangle('Position', [xCoor2by2(sp)-0.5 yCoor2by2(sp)-0.5 2 2],'EdgeColor','r','LineStyle','--'); 
        rectangle('Position', [2.5 2.5+n2 5 5],'EdgeColor','r')
        rectangle('Position', [xCoor3by3(sp)-0.5 yCoor3by3(sp)-0.5+n2 3 3],'EdgeColor','g','LineStyle','-.'); 
        rectangle('Position', [xCoor2by2(sp)-0.5 yCoor2by2(sp)-0.5+n2 2 2],'EdgeColor','r','LineStyle','--'); 
%        title(sprintf('gaus:  I:%.02f, x:%.02f ,y:%.02f, sgx:%.02f,sgy:%.02f',cc1disp(2),cc1disp(3),cc1disp(5),cc1disp(4),cc1disp(6)));
        htt = title('fit to peak')
        set(htt,'Color',[1 1 1]);
        %% FIG 5
        subplot(nFigY,nFigX,5); 
        isProfile = 0;
        if isProfile
            lev = 2;
            plot(Ifit1SHOW(:,pixMaxX(sp),sp),'r-.');hold on; plot(Ifit2SHOW(:,pixMaxX(sp),sp),'r'); plot(spotWin(:,pixMaxX(sp),sp),'r-','LineWidth',2); 
            plot(Ifit1SHOW(pixMaxY(sp),:,sp)+lev,'g-.');hold on; plot(Ifit2SHOW(pixMaxY(sp),:,sp)+lev,'g-'); plot(spotWin(pixMaxY(sp),:,sp)+lev,'g-','LineWidth',2); hold off;
            title(sprintf('cross section X(green) Y(red)',[])); legend('1-GausFitY','2-GausFitY','dataY','1-GausFitX','2-GausFitX','dataX','Location','SouthOutside');
        else
            imgDiff = (spotWin(:,:,sp) - Ifit2SHOW);
            %mnImg = min(imgDiff(:)); mxImg = max(imgDiff(:));
            %imgDiff = (imgDiff - mnImg)/(mxImg-mnImg)*63+1;
            imagesc(imgDiff(y7y1:y7y2,x7x1:x7x2));
            set(gca,'XColor',[1 1 1],'YColor',[1 1 1])
            gcaImg(4) = gca;
            
            axis image
            hC = colorbar; 
            %[ticks64, tickPos, ticks ] = plotColorBar(imgDiff,6,0);
            %set(hC,'Ytick',tickPos,'YTicklabel',ticks);
        end
        htt = title('difference between image and fit : image - fit');
        set(htt,'Color',[1 1 1]);
        %% FIG 6
        subplot(nFigY,nFigX,6); 
        binnedFrm = frm(sp)-1;
        clear temp
        for iAv = 1 : binFrame
            frmRead = binnedFrm*binFrame+iAv;
            temp(:,:,iAv) = double(imread(fname,frmRead)); 
        end        
        img1 = mean(temp(yy1:yy2,xx1:xx2,:),3);
        hImg = imagesc(img1); % axis image; 
        gcaImg(5) = gca;
        axis image; hold on; scatter(xBW(sp),yBW(sp),'w'); hold off
        htt = title('whole image');
        set(htt,'Color',[1 1 1]);

        %% FIG 7-8: neighbouring frames
        subplot(nFigY,nFigX,[7 8]); 
        if frm(sp) > 1
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
        gcaImg(6) = gca;
        
        [szY,szX]=size(IMG);
        
        if frm(sp) > 1
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
        htt = title('whole image with neighbouring frames');
        set(htt,'Color',[1 1 1]);
        
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
    
    % find color boundaries for colorbar normalization
    for i = [1:2,5:6]
        CLim(:,i) = get(gcaImg(i),'CLim');
    end
    CLimMax = max(CLim(2,:));
    CLimMin = min(CLim(1,:));
    CLim = [CLimMin CLimMax];
    for i = [1:2,5:6]
        set(gcaImg(i),'CLim',CLim);
    end
    clear CLim
    for i = [3:4]
        CLim(:,i) = get(gcaImg(i),'CLim');
    end
    CLimMax = max(CLim(2,:));
    CLimMin = min(CLim(1,:));
    CLim = [CLimMin CLimMax];
    for i = [3:4]
        set(gcaImg(i),'CLim',CLim);
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
            [valMin, ixMin] = min(abs(ixDisp - x_ ));
            sp = ixDisp(ixMin);
        end


        while y_ < 1 && y_ > 0 && (gca == hAxInt(1) || gca == hAxInt(2))     % intensity                 
            title(hAx_spotSel1(1),'spots below the intensity threshold, select error threshold to further filter')
            title(hAx_spotSel2(1),'spots below the intensity threshold, select error threshold to further filter')            
            set(hTextAct,'String','select intensity threshold')
            err5 = errFit7by7plot;
            ixInt5 = find(intFitNorm>=y_); % selected
            ixInt5_2 = find(intFitNorm<y_); % some will be selected by Err performance
            set(hLineErr5,      'YData',err5(ixInt5_2));
            set(hLineErr5,      'XData',ixInt5_2);
            set(hLineIntPeak,   'YData',intFitNorm(ixInt5_2));
            set(hLineIntPeak,   'XData',ixInt5_2);
            set(hLineErr5_2,      'YData',err5(ixInt5_2));
            set(hLineErr5_2,      'XData',ixInt5_2);
            set(hLineIntPeak_2,   'YData',intFitNorm(ixInt5_2));
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
            err5 = errFit7by7plot(ixInt5_2);
            ixErr5 = find(err5<=y_); % also selected
            ixErr5_2 = find(err5>y_); % not selected
            ixSel = sort([ixInt5_2(ixErr5) ixInt5]);
            ixNonSel = ixInt5_2(ixErr5_2);
            %ixNonSel = sort([ixErr5_2 ixInt5_2]);

            set(hLineErr5,      'YData',errFit7by7plot(ixSel));
            set(hLineErr5,      'XData',ixSel);
            set(hLineIntPeak,   'YData',intFitNorm(ixSel));
            set(hLineIntPeak,   'XData',ixSel);
            set(hLineErr5_2,      'YData',errFit7by7plot(ixNonSel));
            set(hLineErr5_2,      'XData',ixNonSel);
            set(hLineIntPeak_2,   'YData',intFitNorm(ixNonSel));
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
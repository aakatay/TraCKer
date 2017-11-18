
% 'style','text','BackgroundColor',[1 1 1],
    clear; close all;
    F = findall(0,'type','figure','tag','TMWWaitbar'); delete(F)
    isTraceOnly = 1;
    isCropXY = 0;
    isDebugFilt = 0; 
    isRangeLow = 0;
    noTile = 1;
    binFrameTime = 0; % [ms]
    isCombTraces = 0;
    maxTraceSpeed = nan; % [um/sec] used for calculation of allowed trace jump(not used)
%    minXYspread = 15; % [px] std of the xy location along the trace
    minXYspread = 0.8; % [px] std of the xy location along the trace
    minXYspread = nan;
    minXYspread = inf;
    
%minXYspread = 1.6; % [px] std of the xy location along the trace
    
    
    binFrame = 1;
    
    sptJmpForTracing = 1; % [px]
    %sptReAppearTime = 1; % no gaps allowed
    sptReAppearTime = 2; 
    minTraceLength = 2;
    traceJmpForCombination = 1;   
    
    
    inFocus = 1; % for detection gaussian size. if sharp spots:1, ow:0
    inFocus_w_1X = 0; % if 1.5x mag not used 
    isDiffData = 0;
    
    isBALM = 0;
    isOptimFit = 0;
    isGausFit = 0; % Gaussian localization
    %isGausFitLocalize = 1; % gaussian fit filtering and localization
    StackNum = 1;   
    
    if inFocus_w_1X | isDiffData | inFocus
        WindowSize = 3; 
        BigWindowSize=WindowSize+2;
    else
        WindowSize = 5; 
        %BigWindowSize=WindowSize+4;
        BigWindowSize=WindowSize;
    end
    
    coeffFitParamCoverage = 1000; % number of frames for a fit param
%coeffFitParamCoverage = 50; % number of frames for a fit param
    
        
    %% select frames
    frstFrm = 1; % cropping
    frstFrm2 = 1; 
    lastFrm = [];
    frm1 = 1;
    
    %% read info files
    % read cell info
    
    %% READ the ORIGINAL file
    if exist('fname.mat')
        load fname
        if ~exist('fname')
            fname = sprintf('../\%s',fname);
        end
    end
    
    %defineMask;
    
    fnameBck = [fname(1:end-4) '_background.tif']; % default background file
    label = fname(end-6:end);
    if ~exist(fnameBck)
        fTif1 = dir('../*.tif');
        fTif2 = dir('../../*.tif');
        fTif = [fTif1 ;fTif2];
        for i = 1:length(fTif)
            fTif(i).name;
            if ~isempty(strfind(fTif(i).name,label)) && isempty(strfind(fTif(i).name,'MAX'))
                fnameBck = ['../' fTif(i).name];
                if ~exist(fnameBck)
                    fnameBck = ['../../' fTif(i).name];
                end
            end
        end
        if ~exist(fnameBck)
            fnameBck = fname;
        end
    else % no cropping 
        imgTemp = imread(fnameBck,1);
        BWselect = ones(size(imgTemp));
        save('BWselect','BWselect');
    end
    
    posData_File =dir('posData-coeff*');
    
    %% data filename and image info
    if ~exist('fname')
        ch1_1=strfind(fname,'\'); ch1_2=strfind(fname,'/');
        if length(ch1_1) > length(ch1_2), ch1 = ch1_1; else ch1 = ch1_2; end;
        ch2=strfind(fname,'x');
        ch3=strfind(fname,'xy');
        En1= str2num(fname(ch1(end)+1:ch2(1)-1));
        Boy1= str2num(fname(ch2(1)+1:ch3(1)-1));
    elseif ~isempty(posData_File)
        imgFrst = imread(fname);
        [Boy1,En1]=size(imgFrst);  
        load(posData_File.name,'Frames');
    else
        imgFrst = imread(fname);
        [Boy1,En1]=size(imgFrst);    
        dirFN = dir(fname);
        if dirFN.bytes > 1e10
            Frames=findNumFrames(fname);
        else
            infFN = imfinfo(fname);
            Frames = numel(infFN);
        end
        clear infFN;
            lastFrm = Frames;
        
        Frames = floor(Frames);
    end
    
    if isCropXY
        xx1 = 30; xx2 = 65; % crop
        yy1 = 20; yy2 = 65; % crop
        En1 = xx2-xx1+1;
        Boy1 = yy2-yy1+1;
    else
        xx1 = 1; xx2 = En1; % crop
        yy1 = 1; yy2 = Boy1; % crop
    end    
    sizeImg=[Boy1,En1];
    
    % crop number
    CD=cd;
    if isunix,  slashPos= strfind(CD,'/'); else slashPos= strfind(CD,'\'); end;
    folderName=CD(slashPos(end)+1:end);
    cropNo=sscanf(folderName,'crop%i');
    
    
    coeffMat = 'CoeffFit.mat';
        load(coeffMat);
    Coeff = 0;
    
    assignFileNames


    
    
    
    
    
    
    


    
    %% PLOT2 : recruitment movie
    %% load pre-image
    isTraceOnly = 1;

    %% select ROI
    isCropXY = 0;
    if isCropXY
        preImg = imread(imgFile);
        if exist('imROI.mat')
            load imROI;
            %xx1 = 30; xx2 = 65; % crop
            %yy1 = 20; yy2 = 65; % crop
        else % select the ROI
            imagesc(preImg);
            hPoly = imrect;
            pos = getPosition(hPoly);
            pos = round(pos);
            posRectCrop = pos;
            xx1 = pos(1);
            yy1 = pos(2);
            xx2 = pos(1)+pos(3)-1; 
            yy2 = pos(2)+pos(4)-1;
            save('imROI','xx1','xx2','yy1','yy2');
            %delete(hPoly)
        end
        En1 = xx2-xx1+1;
        Boy1 = yy2-yy1+1;
    else
        xx1 = 1; xx2 = En1; % crop
        yy1 = 1; yy2 = Boy1; % crop
    end    
    
    %% crop trace data 
    load(xyzDataGausFileNm); 
    x=X; y=Y; f=frmNoSpot;
    load(traceDataFileNm0)
    load(traceJmplessDataFileNm)
    save(traceDataFileNm0,'TraceX','TraceY','TraceINT','trInf','frmNoTrace'); % update cfg
        
    TraceX2 = TraceX2 - xx1 + 1;
    TraceY2 = TraceY2 - yy1 + 1;
    trInf(:,4) = trInf(:,4) - xx1+1;
    trInf(:,5) = trInf(:,5) - yy1+1;
    ROItracesX = (trInf(:,4)>1) .* (trInf(:,4)<=En1);
    ROItracesY = (trInf(:,5)>1) .* (trInf(:,5)<=Boy1); 
    ROItraces = ROItracesX .* ROItracesY;
    trInf = trInf(find(ROItraces),:); % select traces in ROI
    
    
    
    %% filter out short traces
    %minTrLenDisp = 2;
    %if minTrLenDisp > 2 
    trInf3 = trInf(trInf(:,2)>=3,:); % select long traces
    trInf4 = trInf(trInf(:,2)>=4,:); % select long traces
    %end
        
    %% filter out spread traces (wondering molecules)
    if ~exist('../stats/'),mkdir('../stats/');end
    hist(trInf(:,8),floor(max(trInf(:,8)*10))); % histogram of spreading of the traces
    set(gca,'Units','pixels'); ylim=get(gca,'Ylim');
    %hline = line([minXYspread, minXYspread], [0,ylim(2)]);
    xlabel('trace spreading [px]');
    ylabel('# of traces')
    title(sprintf('trace spreading histogram. filter threshold:%.02f ',minXYspread));
    imgFig = getframe(gcf); imgOut = imgFig.cdata;
    %set(hline,'Color',[1 0 0])
    imwrite(imgOut,'../stats/traceSpreading.tif')    
    trInf = trInf(trInf(:,8)<=minXYspread,:); % discard spread traces
    trInf3 = trInf3(trInf3(:,8)<=minXYspread,:); % discard spread traces
    trInf4 = trInf4(trInf4(:,8)<=minXYspread,:); % discard spread traces
    
    %% generate trace center array
    TraceCx = nan(size(TraceX));
    TraceCy = nan(size(TraceX));
    for i = 1:size(trInf,1) 
        trix = trInf(i,3):trInf(i,3)+trInf(i,2)-1;
        TraceCx(trix) = trInf(i,4);
        TraceCy(trix) = trInf(i,5);
    end
    %% get bleaching data
    posData_File =dir('posData-coeff*');
    load(posData_File.name,'IMGmean','IMGmax','frmImgMax')
    IMGintNorm = max(IMGmean)./IMGmean;
    %maximg = max(IMGmax(frstFrm2:end));
    %[maximg frmImgMax] = max(IMGmax(frstFrm2:end));
[maximg frmImgMax] = max(IMGmax(1:end));
    %frstFrm2 = cfg.img.frstFrm;
    %maximg = mean(IMGmax(frstFrm2:end));%+std(IMGmax(frstFrm2:end));
    maxNorm = maximg.*IMGintNorm(frmImgMax); % normalized maximum intensity
    %maxNorm = maximg;
    IMGintNorm = smooth(IMGintNorm,500);
    
    %% PLOT3 : color coding & padding % image frame
    isCropData = 0; % frames
    if isCropData
        Tr1=1;
        Tr2=300;
        %Tr2=size(TraceY2,1);
        tt1 = 1; 
        tt2 = size(TraceY2,2);
        tt2 = 50;
        dTr = 1;
        TraceX2 = TraceX2(Tr1:dTr:Tr2,tt1:tt2);
        TraceY2 = TraceY2(Tr1:dTr:Tr2,tt1:tt2);
    end

    % crop frames
    Frames = numel(ixSptFrm)-1; % number of frames
    
    % color data
    %load(traceDataFileNm0,'TraceX2');
    CplotVecN = size(TraceX2,1); % # of traces
    Nframe = size(TraceX2,2); % # of frames
    useCData = 1;
    CData = 1:Frames; % CData=repmat(CData,[size(TraceX2,1) 1]);
    %TraceX2 = TraceX2_; clear TraceX2_;
    if useCData 
        Zmin = min(CData(CData~=0));
        Zmax = max(CData(:));
        Zrange = Zmax - Zmin;
        Cplot = round(63*(CData-Zmin) / Zrange)+1;
        CM = colormap('jet');
    else
        CplotVec = 1: CplotVecN;
        CplotVec = rem(CplotVec,64)+1;
        Cplot = repmat(CplotVec',1,Nframe);
        CM = colormap('lines');
    end
        
    isCropTr=0;
    if isCropTr
        Ntr0 = 1; Ntr = size(TraceX2,1);
        trVec = [Ntr0:Ntr];
        TraceX2 = TraceX2(trVec,1:Frames);
        TraceY2 = TraceY2(trVec,1:Frames);
    end
    padSize = 0;
    zeroIx = find(X == 0);
    X = X + padSize; Y = Y + padSize;
    X(zeroIx) = 0; Y(zeroIx) = 0;   
    zeroIx = find(TraceX2 == 0);
    TraceX2 = TraceX2 + padSize; TraceY2 = TraceY2 + padSize;
    TraceX2(zeroIx) = 0; TraceY2(zeroIx) = 0;
    [Boy2]=size(TraceX2,1);
    
    nonZeroIx = find(TraceY2>0);
    dispPx = {'*',2}; % plot3k
    minimum = 1;
    maximum = max(Cplot(:));

    % image
    load fname;
    imgZFout = ['TraceImage-' fname(4:end)];
    img2D = imread(fname,1); 
    img2D = img2D(yy1:yy2,xx1:xx2);
    img2D = padarray(img2D,[padSize padSize]);
    imgFrm = uint16(zeros(size(img2D))); % image padding frame
    imgFrm(yy1+3:yy2-3,xx1+3:xx2-3)=1;    
    imSz = size(img2D');
    TraceY2(nonZeroIx) = imSz(2)-TraceY2(nonZeroIx)+1;

    mag = 4; % display size
    [mag, pos, m, n ] = calcMaxMag(img2D,mag);
    imSzMag = [m,n]; % magnified image size
    colormap('gray');
    pos(1) = pos(1) - 800;
    pos(1) = pos(1)- 700;
    distFig = 220;
    if distFig < m*2, distFig = m*2; end;
    pos2 = pos; pos2(1) = pos2(1)+ distFig + 20;
    pos3 = pos2;
    pos3(1) = pos3(1)+ distFig + 20;
    pos4 = pos3;
    pos4(1) = pos4(1)+ distFig + 20;
    pos5 = pos4;
    pos5(1) = pos5(1)+ distFig + 20;

    %% PLOT4 : draw TRACES
    NAN = find(isnan(TraceX2));
    TraceX2(NAN)=0;
    TraceY2(NAN)=0;
    
    % find the frames where the traces disappear
    
    hat = ones(1,sptReAppearTime); 
    dspTrcFrm = zeros(size(trInf,1),1);
    for tr = 1:size(trInf,1) % all traces
        dspTrcFrm(tr) = trInf(tr,1)+trInf(tr,2)-1; % last frame
    end
    
    

    hQ = 0; %hImg = image; 
%     lastX = TraceX2(trInf(trInf(:,1)==1,3)); % x position of traces in the first frame
%     lastY = TraceY2(trInf(trInf(:,1)==1,3));
    ix1 = find(trInf(:,1)==1); 
    frm1 = 1;
    tit = 'image';
    m_ = uint16(imSz(1)); n_ = uint16(imSz(2));
    %pos=get(0,'ScreenSize');
    %pos=uint16(pos(3:4)) - [m n-35];
    figImg = figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos2/2 m n]);
    axeImg = axes('Parent',figImg,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren','XLim',0.5+[0 m],'YLim',0.5+[0 n]);
    if ~isTraceOnly
        figImg2 = figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos3/2 m n]);
        axeImg2 = axes('Parent',figImg2,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren','XLim',0.5+[0 m],'YLim',0.5+[0 n]);
        figImg3 = figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos4/2 m n]);
        axeImg3 = axes('Parent',figImg3,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren','XLim',0.5+[0 m],'YLim',0.5+[0 n]);
        figImg4 = figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos5/2 m n]);
        axeImg4 = axes('Parent',figImg4,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren','XLim',0.5+[0 m],'YLim',0.5+[0 n]);
    end
    hWB =  waitbar(0,'marking spots...');
    fig=figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos/2 m n]);
    axe=axes('Parent',fig,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren','XLim',0.5+[0 m],'YLim',0.5+[0 n]);
    ixOld = zeros(size(trInf,1),1);
    %frameVec = 1:6000;
    % image generation parameters
    pxMag = 4; % pixel size scale (from recruitmentTrack)
    imSzBin = imSz*pxMag;
    binImg = zeros([imSzBin(2) imSzBin(1)]); % high res image
    binImgTrSum = binImg; 
    binImgTrSum3 = binImg; 
    binImgTrSum4 = binImg; 
    
    % trace image
    fr = trInf(:,1);
    fr3 = trInf3(:,1);
    fr4 = trInf4(:,1);
    fr(fr==1)=2;
    fr3(fr3==1)=2;
    fr4(fr4==1)=2;
    CLIM = [1 maxNorm];
    
    frameVec = 1:Frames-frstFrm2+1;
    frameVec = 1:6000;
    %frameVec = 405:412;
    
    % frame loop : generate TraceImage =======================
    for ixFrm = frameVec(1:end-1)+1 % all frames % starts with 2nd frame
        frmRead = ixFrm+frstFrm2-1; 
        fname2 = sprintf('../%s',fname); % other channel
        img2D = imread(fname,frmRead);
%img2D = imread(fname2,frmRead);
        img2D = flipud(img2D(yy1:yy2,xx1:xx2));
        img2D = IMGintNorm(ixFrm)*img2D; % intensity normalization
        % padding
        img2D = padarray(img2D,[padSize padSize]);
        %img2D = img2D.*imgFrm;
        %axe=axes('Parent',fig,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','XLim',0.5+[0 m],'YLim',0.5+[0 n]);
        %hold on
        figure(fig);
        set(axe,'Parent',fig,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','XLim',1+[0 m_],'YLim',1+[0 n_]);
        
        hImg = imagesc(img2D,'Parent',axe,CLIM); %axis image; 
        axis tight
        
        %% trace data points
        % current traces
        ixCurr = (ixFrm >= trInf(:,1)+1) .* (ixFrm < (trInf(:,1)+trInf(:,2))); % excluding last data point
        ix = find(ixCurr); % current traces
        ixPos = trInf(ix,3)+ixFrm-trInf(ix,1)-1; % index for x-y positions
        % active traces (including end points)
        ixAct = (ixFrm >= trInf(:,1)) .* (ixFrm < (trInf(:,1)+trInf(:,2))); 
        ix_ = find(ixAct); % current traces
        ixPosAct = trInf(ix_,3)+ixFrm-trInf(ix_,1); % index for x-y positions
        

        %% traces
        X = TraceX(ixPosAct);
        Y = TraceY(ixPosAct); 
        gapPos = find(X.*Y == 0);
        X = X - xx1 + 1;
        Y = Y - yy1 + 1;
        X(gapPos) = [];
        Y(gapPos) = [];
        
        X3 = TraceCx(ixPosAct);
        Y3 = TraceCy(ixPosAct); 
        X3 = X3 - xx1 + 1;
        Y3 = Y3 - yy1 + 1;
        X(X<=0)=nan;Y(X<=0)=nan;
        
        %% detections
        
        ixSel = find((ixFrm==f));
        xsel = x(ixSel);
        ysel = y(ixSel);
        xsel0 = xsel;
        ysel0 = ysel;
        isTraceX = ismember(xsel,X);
        isTraceY = ismember(ysel,Y);
        isTrace = find(isTraceX.*isTraceY);
        xsel(isTrace) = [];
        ysel(isTrace) = [];        
        
        %% -1- mark detections
        nextX = TraceX2(ixPos+1);
        nextY = TraceY2(ixPos+1);
        currX = TraceX2(ixPos);
        currY = TraceY2(ixPos);
        uistack(hImg,'bottom');
        hold on
        dspTrcIx = find(~ixCurr.*ixOld); % index for dissappearing traces
        showTrace =1; % puts arrows
        isShowTrace = 0;
        if showTrace && isShowTrace && ixFrm > 1
            hQdel = hQ(dspTrcIx,:); % handles for discont. traces
            delete(hQdel(hQdel~=0));    % remove the traces of the discontinued traces
        end
        isColorbyTime = 1;
        if isShowTrace
            for i = 1:round(length(ix)) % for each trace 
                iL = ix(i);  % index for each line
                hQ(iL,ixFrm-frm1+1)=quiver(currX(i),currY(i),nextX(i)-currX(i),nextY(i)-currY(i),'Color',CM(Cplot(ixFrm),:));
                if ~showTrace
                    adjust_quiver_arrowhead_size(hQ(i,ixFrm-frm1+1),5)
                end
            end
        end
        
        hscat0 = scatter(X,Boy1-Y+1,1,'g','s','LineWidth',1);
        hscat = scatter(xsel,Boy1-ysel+1,1,'c','s','LineWidth',1);
        %hq = quiver(X,(Boy1-Y+1),0.1,0,'g');
        
        if ~isempty(gapPos)
            X2 = TraceX2(ixPosAct(gapPos));
            Y2 = TraceY2(ixPosAct(gapPos));
            hscat2 = scatter(X2,Y2,1,'r','s','LineWidth',1);
        else
            hscat2 = [];
        end
        %hscat3 = scatter(X3,Boy1-Y3+1,1,'m','.','LineWidth',1);
        
        %% -2- generate image (detections)
        pxX = round(X*pxMag);
        pxY = round(Y*pxMag);
        for i = 1:numel(pxY)
            if pxY(i)<1 || pxX(i)<1 || pxY(i)> imSzBin(2) || pxX(i) > imSzBin(1), continue; end;
            binImg(pxY(i),pxX(i))=1 ...
           +binImg(pxY(i),pxX(i));
        end
        % each detection of the events
        figure(figImg); 
        set(axeImg,'Parent',figImg,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','XLim',1+[0 m_*pxMag],'YLim',1+[0 n_*pxMag]);
        %hImg = imagesc(img2D,'Parent',axeImg); %axis image;
        hold on;hold off;
        imagesc(sum(binImg,3))
        axis equal; axis tight
        hold on;
        dd = 0; %0.05;
        scatter(X*pxMag-dd,Y*pxMag-dd,1,'s','g','MarkerFaceColor','g');
        hold off;

        if ~isTraceOnly
            % -3- generate image (each recruitment) minTrLen:2
            ixSel = find((ixFrm==fr));
            TraceXbin = trInf(ixSel,4);
            TraceYbin = trInf(ixSel,5);
            pxX = round(TraceXbin*pxMag);
            pxY = round(TraceYbin*pxMag);
            binImgTr = zeros([imSzBin(2) imSzBin(1)]);
            for i = 1:numel(pxY)
                if pxY(i)<1 || pxX(i)<1 || pxY(i)> imSzBin(2) || pxX(i) > imSzBin(1), continue; end;
                binImgTr(pxY(i),pxX(i))=1 ...
               +binImgTr(pxY(i),pxX(i));
            end        
            figure(figImg2); 
            set(axeImg2,'Parent',figImg2,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','XLim',1+[0 m_*pxMag],'YLim',1+[0 n_*pxMag]);
            hold on;hold off;
            binImgTrSum = binImgTrSum + binImgTr;
            imagesc(sum(binImgTrSum,3))
            axis equal; axis tight
            hold on;
            scatter(pxX,pxY,1,'s','g','MarkerFaceColor','g');
            hold off;

            % -4- generate image (each recruitment) minTrLen:3
            ixSel = find((ixFrm==fr3));
            TraceXbin = trInf3(ixSel,4);
            TraceYbin = trInf3(ixSel,5);
            pxX = round(TraceXbin*pxMag);
            pxY = round(TraceYbin*pxMag);
            binImgTr3 = zeros([imSzBin(2) imSzBin(1)]);
            for i = 1:numel(pxY)
                if pxY(i)<1 || pxX(i)<1 || pxY(i)> imSzBin(2) || pxX(i) > imSzBin(1), continue; end;
                binImgTr3(pxY(i),pxX(i))=1 ...
               +binImgTr3(pxY(i),pxX(i));
            end        
            figure(figImg3); 
            set(axeImg3,'Parent',figImg3,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','XLim',1+[0 m_*pxMag],'YLim',1+[0 n_*pxMag]);
            hold on;hold off;
            binImgTrSum3 = binImgTrSum3 + binImgTr3;
            imagesc(sum(binImgTrSum3,3))
            axis equal; axis tight
            hold on;
            scatter(pxX,pxY,1,'s','g','MarkerFaceColor','g');
            hold off;

            % -5- generate image (each recruitment) minTrLen:4
            ixSel = find((ixFrm==fr4));
            TraceXbin = trInf4(ixSel,4);
            TraceYbin = trInf4(ixSel,5);
            pxX = round(TraceXbin*pxMag);
            pxY = round(TraceYbin*pxMag);
            binImgTr4 = zeros([imSzBin(2) imSzBin(1)]);
            for i = 1:numel(pxY)
                if pxY(i)<1 || pxX(i)<1 || pxY(i)> imSzBin(2) || pxX(i) > imSzBin(1), continue; end;
                binImgTr4(pxY(i),pxX(i))=1 ...
               +binImgTr4(pxY(i),pxX(i));
            end        
            figure(figImg4); 
            set(axeImg4,'Parent',figImg4,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','XLim',1+[0 m_*pxMag],'YLim',1+[0 n_*pxMag]);
            hold on;hold off;
            binImgTrSum4 = binImgTrSum4 + binImgTr4;
            imagesc(sum(binImgTr4,3))
            axis equal; axis tight
            hold on;
            scatter(pxX,pxY,1,'s','g','MarkerFaceColor','g');
            hold off;
        end

        figure(fig);        
        
        %% print images
        imgFig = getframe(fig);
        imgOut = imgFig.cdata;
        figPos = get(fig,'Position');
        %imgOut = imgOut((figPos(4)-m)/2+1:(figPos(4)+m)/2,(figPos(3)-n)/2+1:(figPos(3)+n)/2,:);
        imgOut1 = imgOut((figPos(4)-n)/2+1:(figPos(4)+n)/2,(figPos(3)-m)/2+1:(figPos(3)+m)/2,:);
        
        
        imgFig = getframe(figImg);
        imgOut = imgFig.cdata;
        figPos = get(figImg,'Position');
        %imgOut = imgOut((figPos(4)-m)/2+1:(figPos(4)+m)/2,(figPos(3)-n)/2+1:(figPos(3)+n)/2,:);
        imgOut2 = imgOut((figPos(4)-n)/2+1:(figPos(4)+n)/2,(figPos(3)-m)/2+1:(figPos(3)+m)/2,:);
        
        if ~isTraceOnly
            imgFig = getframe(figImg2);
            imgOut = imgFig.cdata;
            figPos = get(figImg2,'Position');
            %imgOut = imgOut((figPos(4)-m)/2+1:(figPos(4)+m)/2,(figPos(3)-n)/2+1:(figPos(3)+n)/2,:);
            imgOut3 = imgOut((figPos(4)-n)/2+1:(figPos(4)+n)/2,(figPos(3)-m)/2+1:(figPos(3)+m)/2,:);

            imgFig = getframe(figImg3);
            imgOut = imgFig.cdata;
            figPos = get(figImg3,'Position');
            %imgOut = imgOut((figPos(4)-m)/2+1:(figPos(4)+m)/2,(figPos(3)-n)/2+1:(figPos(3)+n)/2,:);
            imgOut4 = imgOut((figPos(4)-n)/2+1:(figPos(4)+n)/2,(figPos(3)-m)/2+1:(figPos(3)+m)/2,:);

            isFig4 = 1;
            if isFig4
                imgFig = getframe(figImg4);
                imgOut = imgFig.cdata;
                figPos = get(figImg4,'Position');
                %imgOut = imgOut((figPos(4)-m)/2+1:(figPos(4)+m)/2,(figPos(3)-n)/2+1:(figPos(3)+n)/2,:);
                imgOut5 = imgOut((figPos(4)-n)/2+1:(figPos(4)+n)/2,(figPos(3)-m)/2+1:(figPos(3)+m)/2,:);
            end
        end
        
        wp = 5; % width of the padding
        mxV = max(max(max([imgOut1 imgOut2])));
        pad = zeros(size(imgOut1,1),wp,3);
        pad(:,:,1) = mxV/5;
        pad(:,:,2) = mxV/5;
        pad(:,:,3) = mxV/2;
        

        imgOut = [imgOut1 pad imgOut2];
        
        if ixFrm == frameVec(2)
            if exist([imgZFout])
                delete([imgZFout]);
            end
            imwrite(imgOut,imgZFout,'Compression', 'none') 
        else
            imwrite(imgOut,imgZFout,'WriteMode','append','Compression', 'none') 
        end
        delete(hImg)
        delete(hscat0)
        delete(hscat)
        %delete(hq); 
        delete(hscat2)
        %delete(hscat3)

        %waitforbuttonpress
        if ~showTrace
            delete(hQ(hQ(:,ixFrm-1)~=0,ixFrm-1));
        end
        ixOld = ixCurr;
        waitbar(ixFrm/Frames)
    end
    close(hWB);
    hold off; % release the figure 
    
    
    
    isWrite2XLS = 0;
    if isWrite2XLS
        
        load(traceDataFileNm0);
        frameBound = [frameVec(1) frameVec(end)];
        frameBound = [405 412];
        seltrInf = (trInf(:,1)>=frameBound(1)).* (trInf(:,1)<=frameBound(2));
        %seltrInf = (trInf(:,1)<=frameBound(1)).* ((trInf(:,1)+trInf(:,2))>=frameBound(1));
        trInf2 = trInf(find(seltrInf),:);
                
        [TraceXarray TraceYarray] = convertLinear2Array(TraceX,TraceY,trInf);
        [bb nn]=find(sum(TraceXarray,1)>0);
        TraceXarray = TraceXarray(:,frameBound(1):end);
        TraceYarray = TraceYarray(:,frameBound(1):end);
        
        load(traceJmplessDataFileNm)
        [TraceXjumpless TraceYjumpless] = convertLinear2Array(TraceX2,TraceY2,trInf);
        TraceXjumpless = TraceXjumpless(:,frameBound(1):end);
        TraceYjumpless = TraceYjumpless(:,frameBound(1):end);
        
        % write to excel
        fname = 'TraceXjumpless.xls';
        xlswrite(fname,TraceXjumpless)
        fname = 'TraceYjumpless.xls';
        xlswrite(fname,TraceYjumpless)
        
        fname = 'TraceX.xls';
        xlswrite(fname,TraceXarray)
        fname = 'TraceY.xls';
        xlswrite(fname,TraceYarray)

    end
    
        
    return;
    
    %tiff2stack('TraceImage'); % combine frames in a stack
    isTrace = 1;
    save('isTraceFile','isTrace');
    recruitmentTrack;
    isTrace = 0;
    save('isTraceFile','isTrace');    
    recruitmentTrack;
    delete('isTraceFile.mat');    
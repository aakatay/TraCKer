% write images of acquision overlayed with trace & detection points
% RUN in
% (1) bulk data : _00*-coeff**** folder 
% (2) realtime data     : main folder 
%    + after converting trace data by dbgTraces_convertRT.m

% trace points
% green : trace point
% red   : missing trace data
% cyan  : detection (not linked with any trace)

    clear; close all;
    F = findall(0,'type','figure','tag','TMWWaitbar'); delete(F)
        
    %% READ the ORIGINAL file
    isRT = 0;
    isRTrealData = 0; % ow walking average data
    isDispDetect = 0;
    if exist('fname.mat')
        load fname
    elseif exist('fname0.mat')
        load fname0
        fname0_ = fname0;
        load('cfgRT');
        ndigit = cfg.ndigit; % # of digits for sequence number
        digitFormat = sprintf('%%0%1ii',ndigit);
        
        if ~isRTrealData
            fname0 = [cfg.outDIR fname0 'WA_'];
        end   
        fname = [ fname0 num2str(1,digitFormat) '.tif'];
        
        isRT = 1;
    end
        
    %% output file
    
    if isRT
        imgZFout = ['TraceImage-' fname0_ 'RT.tif'];
    else
        imgZFout = ['TraceImage-' fname(4:end)];
    end
    %% data filename and image info
  
    imgFrst = imread(fname);
    [Boy1,En1]=size(imgFrst);    
    if isRT
        fn0 = dir([fname0 '*.tif']);
        Frames = numel(fn0);
    else
        infFN = imfinfo(fname);
        Frames = numel(infFN);
        clear infFN;
    end

    
    %% crop coors
    xx1 = 1; xx2 = En1; % crop
    yy1 = 1; yy2 = Boy1; % crop
    if isRT && isRTrealData
        crp = cfg.crop;
        xx1 = crp(1); xx2 = xx1 + crp(3) - 1; % crop
        yy1 = crp(2); yy2 = yy1 + crp(4) - 1; % crop
        En1 = crp(3);
        Boy1 = crp(4);
    end
    
    
            
            
    
    
    
    
    
    


    
    %% PLOT2 : recruitment movie
    %% load pre-image

    %% crop trace data 
    x=[]; y=[]; f=[];
    if isRT % realtime data
        load('traceDataRT.mat')
        TraceX(isnan(TraceX))=0;
        TraceY(isnan(TraceY))=0;
        TraceX2 = TraceX;
        TraceY2 = TraceY;
    else
        xyzDataGausFile = dir('xyzDataGaus-coeff*.mat');
        traceDataFile = dir('traceData0-coeff*.mat');
        traceJmplessDataFile = dir('traceJmplessData-coeff*.mat');
        
        load(xyzDataGausFile.name); 
        x=X; y=Y; f=frmNoSpot;
        load(traceDataFile.name); % trInf
        load(traceJmplessDataFile.name); % TraceX2
    end
    
            
    %% generate trace center array
    %% get bleaching data
    if isRT % realtime data
        maxNorm = max(imgFrst(:));
        IMGintNorm = ones(Frames,1);
    else 
        posData_File =dir('posData-coeff*');
        load(posData_File.name,'IMGmean','IMGmax','frmImgMax')
        IMGintNorm = max(IMGmean)./IMGmean;
        %maximg = max(IMGmax(frstFrm2:end));
        %[maximg frmImgMax] = max(IMGmax(frstFrm2:end));
        [maximg, frmImgMax] = max(IMGmax(1:end));
        %frstFrm2 = cfg.img.frstFrm;
        %maximg = mean(IMGmax(frstFrm2:end));%+std(IMGmax(frstFrm2:end));
        maxNorm = maximg.*IMGintNorm(frmImgMax); % normalized maximum intensity
        %maxNorm = maximg;
        IMGintNorm = smooth(IMGintNorm,500);
    end
    
    imSz = size(imgFrst);
    
    %% select traces
    %trInf = trInf(1:20,:);
    %trInf = trInf(10,:);
    
    
    % flip
    if 1
        nonZeroIx = find(TraceY>0);
        TraceY(nonZeroIx) = imSz(1)-TraceY(nonZeroIx)+1;
        nonZeroIx = find(TraceY2>0);
        TraceY2(nonZeroIx) = imSz(1)-TraceY2(nonZeroIx)+1;
    end
        
    if 0 
        zeroIx = find(X == 0);
        X = X + padSize; Y = Y + padSize;
        X(zeroIx) = 0; Y(zeroIx) = 0;   
        zeroIx = find(TraceX2 == 0);
        TraceX2 = TraceX2 + padSize; TraceY2 = TraceY2 + padSize;
        TraceX2(zeroIx) = 0; TraceY2(zeroIx) = 0;
    end
    


    mag = 4; % display size
    [pxMag, pos, m, n ] = calcMaxMag(imgFrst,mag); % pixel size scale (from recruitmentTrack)
    imSzMag = [m,n]; % magnified image size
    colormap('gray');
    pos(1) = pos(1) - 800;
    pos(1) = pos(1)- 700;
    distFig = 220;
    if distFig < m*2, distFig = m*2; end
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
    
    tit = 'image';
    m_ = uint16(imSz(1)); n_ = uint16(imSz(2));
    % image2 : sum of traces 
    figImg = figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos2/2 m n]);
    axeImg = axes('Parent',figImg,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','XLim',0.5+[0 m],'YLim',0.5+[0 n]);

    % image1: overlay 
    fig=figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos/2 m n]);
    %axe=axes('Parent',fig,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren','XLim',0.5+[0 m],'YLim',0.5+[0 n]);
    axe=axes('Parent',fig,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','XLim',0.5+[0 m],'YLim',0.5+[0 n]);
    
    ixOld = zeros(size(trInf,1),1);
    % image generation parameters
    imSzBin = imSz*pxMag;
    binImg = zeros([imSzBin(2) imSzBin(1)]); % high res image
    
    % trace image
    fr = trInf(:,1);
    fr(fr==1)=2;
    CLIM = [0 maxNorm];
    
    frameVec = 1:Frames;
    
    hWB =  waitbar(0,'marking spots...');
    % frame loop : generate TraceImage =======================
    %for ixFrm = frameVec(1:end-1)+1 % all frames % starts with 2nd 
    for ixFrm = frameVec(1:end) % all frames % starts with 2nd 
        if isRT
            frmRead = 1;
            fname = [fname0 num2str(ixFrm,digitFormat) '.tif'];
        else
            frmRead = ixFrm;
        end
        img2D = imread(fname,frmRead);
        img2D = img2D(yy1:yy2,xx1:xx2);
        img2D = IMGintNorm(ixFrm)*img2D; % intensity normalization
        % padding
        %img2D = img2D.*imgFrm;
        %axe=axes('Parent',fig,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','XLim',0.5+[0 m],'YLim',0.5+[0 n]);
        %hold on
        figure(fig);
        set(axe,'Parent',fig,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','XLim',0.5+[0 m_],'YLim',0.5+[0 n_]);
        hImg = imagesc(img2D,'Parent',axe,CLIM); axis image;    
        %axis tight      
        box off; 
        
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
        %uistack(hImg,'bottom');
        hold on
        dspTrcIx = find(~ixCurr.*ixOld); % index for dissappearing traces
        
        
        hscat0 = scatter(X,Boy1-Y+1,1,'g','s','LineWidth',1);
        if isDispDetect, hscat = scatter(xsel,Boy1-ysel+1,1,'c','s','LineWidth',1); end
        %hq = quiver(X,(Boy1-Y+1),0.1,0,'g');
        
        if ~isempty(gapPos)
            X2 = TraceX2(ixPosAct(gapPos));
            Y2 = TraceY2(ixPosAct(gapPos));
            if X2.*Y2 == 0, X2 =1.5;Y2=1.5;end
            hscat2 = scatter(X2,Boy1-Y2+1,1,'r','s','LineWidth',1);
            cc=3;
        else
            hscat2 = [];
        end
        %hscat3 = scatter(X3,Boy1-Y3+1,1,'m','.','LineWidth',1);
        
        %% -2- generate image (detections)
        pxX = round((X-0.5)*pxMag+0.5);
        pxY = round((Y-0.5)*pxMag+0.5);
        for i = 1:numel(pxY)
            if pxY(i)<1 || pxX(i)<1 || pxY(i)> imSzBin(2) || pxX(i) > imSzBin(1), continue; end
            binImg(pxY(i),pxX(i))=1 ...
           +binImg(pxY(i),pxX(i));
        end
        if ixFrm>4
            %figure(999)
            %imagesc(binImg)
            cc=4;
        end
        % each detection of the events
        figure(figImg); 
        set(axeImg,'Parent',figImg,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','XLim',1+[0 m_*pxMag],'YLim',1+[0 n_*pxMag]);
        %hImg = imagesc(img2D,'Parent',axeImg); %axis image;
        hold on;hold off;
        imagesc(flipud(sum(binImg,3)))
        box off; 
        axis equal; axis tight
        if 0 
            hold on;
            dd = 0; %0.05;
            scatter(X*pxMag-dd,Y*pxMag-dd,1,'s','g','MarkerFaceColor','g');
            hold off;
        end
  
        
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
        
        
        wp = 5; % width of the padding
        mxV = max(max(max([imgOut1 imgOut2])));
        pad = zeros(size(imgOut1,1),wp,3);
        pad(:,:,1) = mxV/5;
        pad(:,:,2) = mxV/5;
        pad(:,:,3) = mxV/2;
        

        imgOut = [imgOut1 pad imgOut2];
        
        if ixFrm == frameVec(1)
            if exist([imgZFout])
                delete([imgZFout]);
            end
            imwrite(imgOut,imgZFout,'Compression', 'none') 
        else
            imwrite(imgOut,imgZFout,'WriteMode','append','Compression', 'none') 
        end
        delete(hImg)
        delete(hscat0)
        if isDispDetect, delete(hscat); end
        delete(hscat2)

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
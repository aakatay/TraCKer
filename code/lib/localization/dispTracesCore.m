function dispTracesCore(varargin)
% write overlay images : 'TraceImage-....tif'
mxIntensity = 9093;
mxIntensity = 8000;
F = findall(0,'type','figure','tag','TMWWaitbar'); delete(F)
nf = 0;
if nargin>0
    nf = varargin{1}; 
end
isUseCropData = 0;
if nargin>1
    isUseCropData = varargin{2}; 
end
    
    %% load data
    % first frame 
    fo = fopen('frames.txt');
    inp = fgetl(fo);
    sc = strfind(inp,':'); % pos. semicolon
    frm1 = str2num(inp(1:sc-1));
    
    if isUseCropData
        fn2 = rdir('**\traceData0-coeff*.mat');
        fn3 = rdir('**\traceJmplessData-coeff*.mat'); 
        load(fn2(1).name); % TraceX TraceY 
        load(fn3(1).name); % TraceX2 TraceY2 
        load traceDataCrops; % trInf 
        fname = 'acqPitCrops.tif';
        frm1 = 1;
    else
        fn1 = rdir('**\xyzDataGaus-coeff*.mat');
        %fn2 = rdir('**\traceData0-coeff*.mat');
        fn3 = rdir('**\traceJmplessData-coeff*.mat');
        load(fn1(1).name);
        %load(fn2(1).name);
        load('traceData_recTrack');
        load(fn3(1).name);
        load fname;
        trInf(:,13:14) = 0; % no xy offset needed
        x=X; y=Y; f=frmNoSpot;
    end
    
    iminf = imfinfo(fname);
    En1 = iminf.Width;
    Boy1 = iminf.Height;

    
    %% PLOT2 : recruitment movie
    xx1 = 1; xx2 = En1; % crop
    yy1 = 1; yy2 = Boy1; % crop

    
    %% filter out short traces
    trInf3 = trInf(trInf(:,2)>=3,:); % select long traces
    trInf4 = trInf(trInf(:,2)>=4,:); % select long traces
        
    %% filter out spread traces (wondering molecules)
    minXYspread = inf;
    if ~exist('../stats/'),mkdir('../stats/');end
    hist(trInf(:,8),floor(max(trInf(:,8)*10))); % histogram of spreading of the traces
    set(gca,'Units','pixels'); ylim=get(gca,'Ylim');
    hline = line([minXYspread, minXYspread], [0,ylim(2)]);
    xlabel('trace spreading [px]');
    ylabel('# of traces')
    title(sprintf('trace spreading histogram. filter threshold:%.02f ',minXYspread));
    imgFig = getframe(gcf); imgOut = imgFig.cdata;
    set(hline,'Color',[1 0 0])
    imwrite(imgOut,'../stats/traceSpreading.tif')    
    trInf = trInf(trInf(:,8)<=minXYspread,:); % discard spread traces
    trInf3 = trInf3(trInf3(:,8)<=minXYspread,:); % discard spread traces
    trInf4 = trInf4(trInf4(:,8)<=minXYspread,:); % discard spread traces
    
    %% generate trace center array
    TraceCx = nan(size(TraceX));
    TraceCy = nan(size(TraceX));
    for i = 1:size(trInf,1) 
        trix = trInf(i,3):trInf(i,3)+trInf(i,2)-1;
        TraceCx(trix) = trInf(i,9);
        TraceCy(trix) = trInf(i,10);
    end
    
    %% PLOT3 : color coding & padding % image frame

    % crop frames
    if nf ~= 0
        Frames = nf;
    else
        Frames = numel(ixSptFrm)-1; % number of frames
    end
    
    % color data
    CplotVecN = size(TraceX2,1); % # of traces
    Nframe = size(TraceX2,2); % # of frames
    useCData = 1;
    CData = 1:Frames; 
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
        
    nonZeroIx = find(TraceY2>0);
    dispPx = {'*',2}; % plot3k
    minimum = 1;
    maximum = max(Cplot(:));

    % image
    imgZFout = ['TraceImage-' fname(4:end)];
    img2D = imread(fname,1); 
    img2D = img2D(yy1:yy2,xx1:xx2);
    imSz = size(img2D');

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
    sptReAppearTime = cfg.trace.sptReAppearTime; %(frames) use the value from tracker function generating trace values
    % find the frames where the traces disappear
    
    hat = ones(1,sptReAppearTime); 
    dspTrcFrm = zeros(size(trInf,1),1);
    for tr = 1:size(trInf,1) % all traces
        dspTrcFrm(tr) = trInf(tr,1)+trInf(tr,2)-1; % last frame
    end
    
    

    hQ = 0; 
    ix1 = find(trInf(:,1)==1); 
    tit = 'image';
    m_ = uint16(imSz(1)); n_ = uint16(imSz(2));
    figImg = figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos2/2 m n]);
    axeImg = axes('Parent',figImg,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren','XLim',0.5+[0 m],'YLim',0.5+[0 n]);
    figImg2 = figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos3/2 m n]);
    axeImg2 = axes('Parent',figImg2,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren','XLim',0.5+[0 m],'YLim',0.5+[0 n]);
    figImg3 = figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos4/2 m n]);
    axeImg3 = axes('Parent',figImg3,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren','XLim',0.5+[0 m],'YLim',0.5+[0 n]);
    figImg4 = figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos5/2 m n]);
    axeImg4 = axes('Parent',figImg4,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren','XLim',0.5+[0 m],'YLim',0.5+[0 n]);
    hWB =  waitbar(0,'marking spots...');
    fig=figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos/2 m n]);
    axe=axes('Parent',fig,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren','XLim',0.5+[0 m],'YLim',0.5+[0 n]);
    ixOld = zeros(size(trInf,1),1);
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
    
    frameVec = 1:Frames-frm1+1;
    %frameVec = 168:Frames-frm1+1;
    
    % frame loop : generate TraceImage =======================
    isDbgTraceLen = 0;
    for ixFrm = frameVec(1:end-1)+1 % all frames % starts with 2nd frame
        frmRead = ixFrm+frm1-1; 
        fname2 = sprintf('../%s',fname); % other channel
        img2D = imread(fname,frmRead);
        img2D = double(flipud(img2D(yy1:yy2,xx1:xx2)));
        %img2D = IMGintNorm(ixFrm)*img2D; % intensity normalization
        
        
        figure(fig);
        set(axe,'Parent',fig,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','XLim',1+[0 m_],'YLim',1+[0 n_]);
        CM = gray(256);
        colormap(CM);
        img2D = mxIntensity-img2D;
        img2D = img2D/mxIntensity*256;
        max(img2D(:));
        
        img2D = uint16(img2D);
        hImg = image(img2D,'Parent',axe); %axis image; 
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
        
        X = TraceX(ixPosAct) + trInf(ix_,13)';
        Y = TraceY(ixPosAct) + trInf(ix_,14)'; 
        gapPos = find(TraceX(ixPosAct).*TraceY(ixPosAct) == 0);
        X = X - xx1 + 1;
        Y = Y - yy1 + 1;
        
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
        nextX = TraceX2(ixPos+1)  + trInf(ix,13)';
        nextY = TraceY2(ixPos+1)  + trInf(ix,14)';
        currX = TraceX2(ixPos)  + trInf(ix,13)';
        currY = TraceY2(ixPos) + trInf(ix,14)';
        uistack(hImg,'bottom');
        hold on
        dspTrcIx = find(~ixCurr.*ixOld); % index for dissappearing traces
        showTrace =0; % puts arrows
        isShowTrace = 0;
        if showTrace && isShowTrace && ixFrm > 1
            hQdel = hQ(dspTrcIx,:); % handles for discont. traces
            delete(hQdel(hQdel~=0));    % remove the traces of the discontinued traces
        end
        isColorbyTime = 1;
        if isShowTrace
            for i = 1:round(length(ix)) % for each trace 
                iL = ix(i);  % index for each line
                hQ(iL,ixFrm)=quiver(currX(i),currY(i),nextX(i)-currX(i),nextY(i)-currY(i),'Color',CM(Cplot(ixFrm),:));
                if ~showTrace
                    adjust_quiver_arrowhead_size(hQ(i,ixFrm),5)
                end
            end
        end
        hscat0 = scatter(X,Boy1-Y+1,1,'g','s','LineWidth',1);
        hscat = scatter(xsel,Boy1-ysel+1,1,'c','s','LineWidth',1);
        %hq = quiver(X,(Boy1-Y+1),0.1,0,'g');
        
        if ~isempty(gapPos)
            X2 = TraceX2(ixPosAct(gapPos)) + trInf(ix_(gapPos),13)';
            Y2 = TraceY2(ixPosAct(gapPos)) + trInf(ix_(gapPos),14)';
            hscat2 = scatter(X2,Boy1-Y2+1,1,'r','s','LineWidth',1);
        else
            hscat2 = [];
        end
        %hscat3 = scatter(X3,Boy1-Y3+1,1,'m','.','LineWidth',1);
        
        %% -2- generate image (detections)
        if 0 % trace detections ow all detections
            xsel = X;
            ysel = Y;
        end
        pxX = round(xsel0*pxMag);
        pxY = round(ysel0*pxMag);
        for i = 1:numel(pxY)
            if pxY(i)<1 || pxX(i)<1 || pxY(i)> imSzBin(2) || pxX(i) > imSzBin(1), continue; end;
            binImg(pxY(i),pxX(i))=1 ...
           +binImg(pxY(i),pxX(i));
        end
        % each detection of the events
        figure(figImg); 
        set(axeImg,'Parent',figImg,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','XLim',1+[0 m_*pxMag],'YLim',1+[0 n_*pxMag]);
        hold on;hold off;
        img2 = sum(binImg,3);
        img2 = img2/max(img2(:))*256;
        img2(isnan(img2)) = 0;
        image(img2);
        axis equal; axis tight
        hold on;
        dd = 0; %0.05;
        scatter(X*pxMag-dd,Y*pxMag-dd,1,'s','g','MarkerFaceColor','g');
        scatter(xsel*pxMag,ysel*pxMag,1,'c','s','LineWidth',1);
        hold off;
        
        % -3- generate image (each recruitment) minTrLen:2
        ixSel = find((ixFrm==fr));
        TraceXbin = trInf(ixSel,9);
        TraceYbin = trInf(ixSel,10);        
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
        img3 = sum(binImgTrSum,3);
        img3 = img3/max(img3(:))*256;
        img3(isnan(img3)) = 0;
        image(img3);
        axis equal; axis tight
        hold on;
        scatter(pxX,pxY,1,'s','g','MarkerFaceColor','g');
        hold off;
        
        % -4- generate image (each recruitment) minTrLen:3
        ixSel = find((ixFrm==fr3));
        TraceXbin = trInf3(ixSel,9);
        TraceYbin = trInf3(ixSel,10);
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
        TraceXbin = trInf4(ixSel,9);
        TraceYbin = trInf4(ixSel,10);
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

        figure(fig);        
        
        %% print images
        imgFig = getframe(fig);
        imgOut = imgFig.cdata;
        figPos = get(fig,'Position');
        imgOut1 = imgOut((figPos(4)-n)/2+1:(figPos(4)+n)/2,(figPos(3)-m)/2+1:(figPos(3)+m)/2,:);
        
        
        imgFig = getframe(figImg);
        imgOut = imgFig.cdata;
        figPos = get(figImg,'Position');
        imgOut2 = imgOut((figPos(4)-n)/2+1:(figPos(4)+n)/2,(figPos(3)-m)/2+1:(figPos(3)+m)/2,:);
        
        
        imgFig = getframe(figImg2);
        imgOut = imgFig.cdata;
        figPos = get(figImg2,'Position');
        imgOut3 = imgOut((figPos(4)-n)/2+1:(figPos(4)+n)/2,(figPos(3)-m)/2+1:(figPos(3)+m)/2,:);
        
        imgFig = getframe(figImg3);
        imgOut = imgFig.cdata;
        figPos = get(figImg3,'Position');
        imgOut4 = imgOut((figPos(4)-n)/2+1:(figPos(4)+n)/2,(figPos(3)-m)/2+1:(figPos(3)+m)/2,:);
        
        isFig4 = 1;
        if isFig4
            imgFig = getframe(figImg4);
            imgOut = imgFig.cdata;
            figPos = get(figImg4,'Position');
            imgOut5 = imgOut((figPos(4)-n)/2+1:(figPos(4)+n)/2,(figPos(3)-m)/2+1:(figPos(3)+m)/2,:);
        end
        
        wp = 5; % width of the padding
        mxV = max(max(max([imgOut1 imgOut2 imgOut3])));
        pad = zeros(size(imgOut3,1),wp,3);
        pad(:,:,1) = 21;
        pad(:,:,2) = 100;
        pad(:,:,3) = 21;
        
        if isDbgTraceLen
            if isFig4
                imgOut = [imgOut1 pad imgOut2 pad imgOut3 pad imgOut4 pad imgOut5];
            else
                imgOut = [imgOut1 pad imgOut2 pad imgOut3 pad imgOut4];
            end
        else
                imgOut = [imgOut1 pad imgOut3];
                imgOut = [imgOut1 pad imgOut2 pad imgOut3];
        end
        
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
        if exist('hscat3'),delete(hscat3); end;

        %waitforbuttonpress
        if isShowTrace && ~showTrace
            delete(hQ(hQ(:,ixFrm-1)~=0,ixFrm-1));
        end
        ixOld = ixCurr;
        waitbar(ixFrm/Frames)
    end
    close(hWB);
    hold off; % release the figure 
end
% run to mark recruitment traces on live acq movie and recruitment movie
    
    clear all;
    close all;
    
    posFN=rdir('**\*posData-coeff0*.mat');
    load(posFN.name, 'IMGmax','frstFrm2')
    [movieMax,~] = max(IMGmax);
    CLIM = [1 movieMax];
    
    avgFN=rdir('AVG*.tif');
    A = imread(avgFN.name);
    %% load trace data 
    load('traceData_recTrack.mat');
    
            
    %% filter out spread traces (wondering molecules)
    if ~exist('stats/'),mkdir('stats/');end
    hist(trInf(:,8),floor(max(trInf(:,8)*10))); % histogram of spreading of the traces
    set(gca,'Units','pixels'); ylim=get(gca,'Ylim');
    xlabel('trace spreading [px]');
    ylabel('# of traces')
    title(sprintf('trace spreading histogram. filter threshold:nan'));
    imgFig = getframe(gcf); imgOut = imgFig.cdata;
    imwrite(imgOut,'stats/traceSpreading.tif')    
    if 0
        trInf = trInf(trInf(:,8)<=minXYspread,:); % discard spread traces
        trInf3 = trInf3(trInf3(:,8)<=minXYspread,:); % discard spread traces
        trInf4 = trInf4(trInf4(:,8)<=minXYspread,:); % discard spread traces
    end
    
    
    %% PLOT3 : color coding & padding % image frame

    Frames = numel(ixSptFrm)-1; % number of frames
    
    % color data
    %load(traceDataFileNm0,'TraceX');
    CplotVecN = size(TraceX,1); % # of traces
    Nframe = size(TraceX,2); % # of frames
    useCData = 1;
    CData = 1:Frames; % CData=repmat(CData,[size(TraceX,1) 1]);
    %TraceX = TraceX_; clear TraceX_;
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
        
    [Boy2]=size(TraceX,1);
    
    nonZeroIx = find(TraceY>0);
    dispPx = {'*',2}; % plot3k

    % image
    load fname;
    tit = 'image';
    imgZFout = ['TraceImage-' fname(4:end)];
    img2D = imread(fname,1); 
    imgFrm = uint16(zeros(size(img2D))); % image padding frame
    imSz = size(img2D');
    TraceY(nonZeroIx) = imSz(2)-TraceY(nonZeroIx)+1;
    [Boy1, En1] = size(img2D);
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

    %% PLOT4 : draw TRACES
    NAN = find(isnan(TraceX));
    TraceX(NAN)=0;
    TraceY(NAN)=0;


    hQ = 0; %hImg = image; 
%     lastX = TraceX(trInf(trInf(:,1)==1,3)); % x position of traces in the first frame
%     lastY = TraceY(trInf(trInf(:,1)==1,3));
    ix1 = find(trInf(:,1)==1); 
    m_ = uint16(imSz(1)); n_ = uint16(imSz(2));
    %pos=get(0,'ScreenSize');
    %pos=uint16(pos(3:4)) - [m n-35];
    figImg = figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos2/2 m n]);
    axeImg = axes('Parent',figImg,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren','XLim',0.5+[0 m],'YLim',0.5+[0 n]);
    figImg2 = figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos3/2 m n]);
    axeImg2 = axes('Parent',figImg2,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren','XLim',0.5+[0 m],'YLim',0.5+[0 n]);
    figImg3 = figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos4/2 m n]);
    axeImg3 = axes('Parent',figImg2,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren','XLim',0.5+[0 m],'YLim',0.5+[0 n]);

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
    fr2 = trInf(:,2);
    %fr(fr==1)=2;
    
    AR = zeros([size(img2D)*pxMag 3]); % recruitment intensity overlay
    AR(:,:,1) = repelem(A,pxMag,pxMag);
    mxRec=5;
    mxA = double(max(A(:)));
    ARscale= double(mxA)/mxRec;
    
    frameVec = 1:Frames;
    ixPlotted = [];
    % frame loop : generate TraceImage =======================
    for ixFrm = frameVec(1:end) % all frames % starts with 2nd frame

        %% trace coors
        ixSel = find((ixFrm>=fr).*(ixFrm<=fr+fr2-1));
        TraceXbin = trInf(ixSel,9)-0.5;
        TraceYbin = trInf(ixSel,10)-0.5;
        pxX = round(TraceXbin*pxMag);
        pxY = round(TraceYbin*pxMag);
                 
        %% images
        % intensity image 1/2
        figure(fig);
        frmRead = ixFrm+frstFrm2-1; 
        img2D = imread(fname,frmRead);
        %img2D = flipud(img2D);
        img2D = repelem(img2D,pxMag,pxMag);
        hImg = imagesc(img2D,CLIM); %axis image; 
        
        % recruitment image 2/2
        figure(figImg2); 
        binImgTr = zeros([imSzBin(2) imSzBin(1)]);
        for i = 1:numel(pxY)
            if ismember(ixSel(i),ixPlotted), continue; end;
            if pxY(i)<1 || pxX(i)<1 || pxY(i)> imSzBin(2) || pxX(i) > imSzBin(1), continue; end;
            binImgTr(pxY(i),pxX(i))=1 ...
           +binImgTr(pxY(i),pxX(i));
        end     
        binImgTrSum = binImgTrSum + binImgTr;
        imagesc(sum(binImgTrSum,3))
        
        if ~isempty(find(binImgTrSum>mxRec)), warning('image saturated: use higher mxRec value'); end
        % recruitment image 3/3
        figure(figImg3); 
        
        
               
        %% display trace coors
        figure(fig);
        hold on;
        scatter(pxX,pxY,1,'.','g','MarkerFaceColor','g');
        scatter(pxX,pxY,100,'o','r');
        hold off;
        set(gca,'Parent',fig,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','XLim',1+[0 m_],'YLim',1+[0 n_]);
        axis equal; axis tight

        figure(figImg2); 
        hold on;
        scatter(pxX,pxY,1,'.','g','MarkerFaceColor','g');
        hold off;
        set(gca,'Parent',figImg2,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','XLim',1+[0 m_*pxMag],'YLim',1+[0 n_*pxMag]);
        axis equal; axis tight
        
        figure(figImg3); 
        AR(:,:,2) = binImgTrSum*ARscale*mxRec/max(binImgTrSum(:)) * 1.2;
        imagesc(AR/mxA*0.7);
        set(gca,'Parent',figImg3,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','XLim',1+[0 m_*pxMag],'YLim',1+[0 n_*pxMag]);
        axis equal; axis tight
        
        %% print images
        imgFig = getframe(fig);
        imgOut = imgFig.cdata;
        figPos = get(fig,'Position');
        %imgOut = imgOut((figPos(4)-m)/2+1:(figPos(4)+m)/2,(figPos(3)-n)/2+1:(figPos(3)+n)/2,:);
        imgOut1 = imgOut((figPos(4)-n)/2+1:(figPos(4)+n)/2,(figPos(3)-m)/2+1:(figPos(3)+m)/2,:);
        if ixFrm == frameVec(1), imgOut1=flipud(imgOut1); end;
            
        imgFig = getframe(figImg2);
        imgOut = imgFig.cdata;
        figPos = get(figImg2,'Position');
        %imgOut = imgOut((figPos(4)-m)/2+1:(figPos(4)+m)/2,(figPos(3)-n)/2+1:(figPos(3)+n)/2,:);
        imgOut3 = imgOut((figPos(4)-n)/2+1:(figPos(4)+n)/2,(figPos(3)-m)/2+1:(figPos(3)+m)/2,:);
        if ixFrm == frameVec(1), imgOut3=flipud(imgOut3); end;
        
        imgFig = getframe(figImg3);
        imgOut = imgFig.cdata;
        figPos = get(figImg3,'Position');
        %imgOut = imgOut((figPos(4)-m)/2+1:(figPos(4)+m)/2,(figPos(3)-n)/2+1:(figPos(3)+n)/2,:);
        imgOut4 = imgOut((figPos(4)-n)/2+1:(figPos(4)+n)/2,(figPos(3)-m)/2+1:(figPos(3)+m)/2,:);
        
        wp = 1; % width of the padding
        mxV = max(max(max([imgOut1 imgOut3])));
        pad = zeros(size(imgOut3,1),wp,3);
        pad(:,:,1) = mxV/5;
        pad(:,:,2) = mxV/5;
        pad(:,:,3) = mxV/2;
        
        imgOut = [imgOut1 pad imgOut3 pad imgOut4];
        
        if ixFrm == frameVec(1)
            if exist([imgZFout])
                delete([imgZFout]);
            end
            imwrite(imgOut,imgZFout,'Compression', 'none') 
        else
            imwrite(imgOut,imgZFout,'WriteMode','append','Compression', 'none') 
        end
        %delete(hImg)

        %waitforbuttonpress
        
        ixPlotted = [ixPlotted ixSel];
        waitbar(ixFrm/Frames)
    end
    close(hWB);
    
    
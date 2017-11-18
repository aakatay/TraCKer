% combines the dorsal and ventral data
% run in Rotated3D folder
close all

Fdorsal = 'dorsal3D\-coeff57';
Fventral = 'ventral3D\-coeff40';
load([Fventral '\traceData-coeff40.mat']);
    TraceXventral = TraceX;
    TraceYventral = TraceY;
    TraceZventral = TraceZ;
    TraceSpeedventral = TraceSpeed;
    TraceINTventral = TraceINT;
load([Fdorsal '\traceData-coeff57.mat']);
    TraceZ = TraceZ+29;
    
    TraceX = [TraceX ; TraceXventral];
    TraceY = [TraceY ; TraceYventral];
    TraceZ = [TraceZ ; TraceZventral];
    TraceSpeed = [TraceSpeed ; TraceSpeedventral];
    TraceINT = [TraceINT ; TraceINTventral];

    load([Fdorsal '\posData-coeff57.mat'],'En1'); %(frames) use the value from tracker function generating trace values
    En1dorsal = En1;
    load([Fventral '\posData-coeff40.mat'],'En1'); %(frames) use the value from tracker function generating trace values
    En1ventral = En1;
    load([Fdorsal '\posData-coeff57.mat'],'Boy1'); %(frames) use the value from tracker function generating trace values
    imSz = [max(En1dorsal,En1ventral),Boy1];
    load([Fdorsal '\posData-coeff57.mat'],'sptReAppearTime'); 
    load([Fdorsal '\inputInfo_wesScope.mat'],'acqRate_sec');
    run2plotDorsalVentral = 1;

save('traceData-combined.mat','TraceX','TraceY','TraceZ','TraceINT','TraceSpeed')
    
    isFewTraces = 0;
    if isFewTraces
        N = 1000;
        TraceX = TraceX(1:N,:);
        TraceY = TraceY(1:N,:);
        TraceZ = TraceZ(1:N,:);
        TraceSpeed = TraceSpeed(1:N,:);
    end
%% CALL the tracker to plot traces
%TraCKer_3D_w_ZcolorPlot_deep_strain
    speedPx2nm = pxSz/acqRate_sec*1000;
    pxSzZ = 0.1/6; % um (isotropic)
    %speedPx2nm=1;
    isHist = 0;
    if isHist
        figure
        hist(TraceSpeed((TraceSpeed(:)>0) )*speedPx2nm,1000);
        title('speed distribution of each trace step')
        xlabel('trace speed (nm/sec)')
        colormap('gray');
    end
    
    % Tracespeed mean
    if 0
        for i = 1:size(TraceSpeed,1) % trace mean
            TraceSpeed(i,TraceSpeed(i,:)>0)= mean(TraceSpeed(i,TraceSpeed(i,:)>0));
        end
    end
    
    if 1
    nAv=3;
        TraceSpeed = full(TraceSpeed);
        for i = 1:size(TraceSpeed,1) % walking average
            TraceSpeed(i,TraceSpeed(i,:)>0)= conv(TraceSpeed(i,TraceSpeed(i,:)>0),ones(nAv,1),'same')/nAv;
        end
    end
    if isHist
        figure
        hist(TraceSpeed((TraceSpeed(:)>0) )*speedPx2nm,1000);
        title('speed distribution of each trace step averaged by neighbours (WAby3)')
        xlabel('trace speed (nm/sec)')
        colormap('gray');
    end
        
    isTestColorCoding = 0;
    if isTestColorCoding
        TraceSpeed = 10:64;
        TraceX = TraceSpeed*3;
        TraceY = TraceSpeed*3;
        TraceSpeed(numel(TraceSpeed)+1)=-0.15;
    end
    
    zMaxHigh = 3.1/pxSzZ;
    %TraceZ(TraceZ>zMaxHigh)=zMaxHigh;
    % color data
    CplotVecN = size(TraceX,1); % # of traces
    Nframe = size(TraceX,2); % # of frames
    CMinZ = 1; % oclor map in Z
    CMinSpeed =0; % color in speed
    %speedMax = 3.5;
    if CMinZ 
        Zmin = min(TraceZ(TraceZ~=0));
        Zmax = max(TraceZ(:));
        Zrange = Zmax - Zmin;
        Cplot = round(63*(TraceZ-Zmin) / Zrange)+1;
        CM = colormap('jet');
        ColorMode = 'Z';
    elseif CMinSpeed
        Zmin = min(TraceSpeed(TraceSpeed~=0));
        Zmax = max(TraceSpeed(:));
        Zrange = Zmax - Zmin;
        Cplot = round(63*(TraceSpeed-Zmin) / Zrange)+1;
        CM = colormap('jet');
        ColorMode = 'Speed';
    else % color by trace
        CplotVec = 1: CplotVecN;
        CplotVec = rem(CplotVec,64)+1;
        Cplot = repmat(CplotVec',1,Nframe);
        CM = colormap('lines');
    end
    Frames = size(TraceX,2);
    En = Frames;
    [Boy2]=size(TraceX,1);
    
    nonZeroIx = find(TraceY>0);
    dispPx = {'*',2}; % plot3k
    minimum = 1;
    maximum = max(Cplot(:));
    
% PLOT
    NAN = find(isnan(TraceX));
    TraceX(NAN)=0;
    TraceY(NAN)=0;
    TraceZ(NAN)=0;

    % parameters:

   
    imgZFout = 'TraceImage';
    magImg = 10;
    mag = 2;

    
    TraceY = imSz(2)-TraceY+1;
    TraceX = TraceX*mag;
    TraceY = TraceY*mag;
    
    figSz(1) = imSz(1)*magImg;
    figSz(2) = imSz(2)*magImg;
    colormap('gray');

    zMax = max(TraceZ(:));
    if zMax == 0, zMax = 1; TraceZ = TraceZ+1; end;

    % find the frames where the traces disappear
    hat = ones(1,sptReAppearTime); 
    dspTrcFrm = zeros(size(TraceX,1),1);
    for i = 1:size(TraceX,1) % all traces
        trace = full(TraceX(i,:));
        trFrm = conv(trace,hat,'same'); % active frames of the trace
        trFrm(1:find(trFrm>0,1)) = 1;
        lastFrm = find(trFrm==0,1); 
        if ~isempty(lastFrm)
            dspTrcFrm(i) = lastFrm;
        end
    end

    imgZFout = 'TraceImageStack';
    imgZFoutScat = 'TraceImage0';
    fOut = sprintf('%s_%s.tif',imgZFout,ColorMode);
    fOutScat = sprintf('%s_%s.tif',imgZFoutScat,ColorMode);
    hQ = 0; %hImg = image; 
    lastX = TraceX(:,1); 
    lastY = TraceY(:,1);
    lastZ = TraceZ(:,1);
    tit = 'image';
    m = imSz(1)*mag; n = imSz(2)*mag;
    pos=get(0,'ScreenSize');
    pos=pos(3:4) - [m n-35];
    fig=figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos/2 m n]);
    axe=axes('Parent',fig,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren','XLim',0.5+[0 m],'YLim',0.5+[0 n]);

    isScatterPlot = 1;
    colormap('jet')
    if isScatterPlot
        xx = TraceX(TraceSpeed>0);
        yy = TraceY(TraceSpeed>0);
        cc = CM(Cplot(TraceSpeed>0),:);
        scatter(xx,yy,4,cc,'*');
        gry = 0.5;
        set(gcf,'Color',[1 1 1]*gry);
        imgFig = getframe(gcf);
        imgOut = imgFig.cdata;
        figPos = get(gcf,'Position');
        %imgOut = imgOut((figPos(4)-m)/2+1:(figPos(4)+m)/2,(figPos(3)-n)/2+1:(figPos(3)+n)/2,:);
        imgOut = imgOut((figPos(4)-n)/2+1:(figPos(4)+n)/2,(figPos(3)-m)/2+1:(figPos(3)+m)/2,:);
        imwrite(imgOut,fOutScat,'Compression', 'none') 
        colorbarFout = sprintf('%s_colorbar.tif',fOutScat(1:end-4));
    else
        colorbarFout = sprintf('%s_colorbar.tif',fOut(1:end-4));

        for ixFrm = 1:En-1
            img2D = imread(fname,ixFrm+1); 
            img2D = flipud(img2D(yy1:yy2,xx1:xx2));
            %delete(hImg);

            hImg = imagesc(img2D,'Parent',axe); %axis image; 
            %set(gca,'position',[0 0 1 1]); set(gcf,'position',[0 0 figSz(1) figSz(2)]);
            ixTrc = find(TraceX(:,ixFrm+1)>0);
            currX = TraceX(:,ixFrm+1);
            currY = TraceY(:,ixFrm+1);
            currZ = TraceZ(:,ixFrm+1);
            uistack(hImg,'bottom'); hold on
            ixTrc2 = find(lastX(ixTrc)>0);
            ix = ixTrc(ixTrc2); % indices of traces tracked in the current frame    
            dspTrcIx = find(ixFrm+1 == dspTrcFrm); % index for dissappearing traces
            showTrace =1;
            if showTrace
                hQdel = hQ(dspTrcIx,:); % handles for discont. traces
                delete(hQdel(find(hQdel~=0)));    % remove the traces of the discontinued traces
            end

            for i = 1:round(length(ix)) % for each trace 
                iL = ix(i);  % index for each line
                hQ(iL,ixFrm)=quiver(lastX(iL),lastY(iL),currX(iL)-lastX(iL),currY(iL)-lastY(iL),'Color',CM(Cplot(iL,ixFrm),:));
                if ~showTrace
                    adjust_quiver_arrowhead_size(hQ(iL,ixFrm-1),5)
                end
            end
            % print images
            imgFig = getframe(gcf);
            imgOut = imgFig.cdata;
            figPos = get(gcf,'Position');
            %imgOut = imgOut((figPos(4)-m)/2+1:(figPos(4)+m)/2,(figPos(3)-n)/2+1:(figPos(3)+n)/2,:);
            imgOut = imgOut((figPos(4)-n)/2+1:(figPos(4)+n)/2,(figPos(3)-m)/2+1:(figPos(3)+m)/2,:);
            if ixFrm == 1
                if exist([imgZFout '.tif'])
                    delete([imgZFout '.tif']);
                end
                imwrite(imgOut,fOut,'Compression', 'none') 
            else
                imwrite(imgOut,fOut,'WriteMode','append','Compression', 'none') 
            end

            delete(hImg)
            %export_fig(imgZFout,'-append');
            fprintf('frame %i/%i \n',ixFrm,En);
            lastX = currX; 
            lastY = currY;
            lastZ = currZ;
        end
        
    end
    
    
    if CMinZ 
        plotColorBar(TraceZ*pxSzZ,5,'Z [\mum]'); % 
    elseif CMinSpeed
        plotColorBar(TraceSpeed*speedPx2nm,5,'speed [nm/sec]');
    end
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'Color',[1 1 1]);
    imgFig = getframe(gcf);
    imgOut = imgFig.cdata;
    imwrite(imgOut,colorbarFout,'Compression', 'none') 
    

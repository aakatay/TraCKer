function traceProj
% Displays the traces of spots on projection movie
% loads the trace info from the 'traceData.mat'
    clear; close all;
    %% data parameters
    tAcq = 2.182; % [sec] cell 1-4
    tAcq = 2.026; % [sec] mu-2
    %tAcq = 2010; % [sec] cell1
    pxSz = 104; %[nm] Wes Scope
    pxSz = 114; %[nm] Bi-Chang Scope
    pxSzZ= 0.4; %[um]
    %pxSzZ= 0.240; %um
    %% display param
    isZinColor = 1; % ow speed in color
    isZinColor3D = 0;
    dispPx = {'o',2}; % plot3k
    CM = ('jet'); %('jet');
    isSpeedinLog = 0;
    upLev = 0; % no speed color leveling
    upLev = 3*pxSz/tAcq;
    lowLev = 0; 
    lowLev = 0/8;
    vibrCorr = 0; % vibr. corrected
    smthWin = 5; % smoothing
    smoothing = 1;
    minLenTrace = 30;
    TraceX = [];
    TraceY = [];
    TraceZ = [];
    TraceINT = [];
    TraceT0 = [];
    

    [traceDataPath, imgFin, imgSpeedFout, imgZFout] = getAddr;
    %traceDataPath = 'traceData.mat'; % vibration corrected

    load(char(traceDataPath)); % load traceData: 'TraceX','TraceY','TraceZ','TraceINT','TraceT0'
    hFig=figure;
    R1 = 1; R2 = 10;
    for iTrace = 1:100
        R1 = iTrace; R2 = iTrace+1;
        plotTraces3D
        pause
        close all
    end
    saveFigImg('outImage.tif',999,hFig)
    return;


    function [traceDataPath, imgFin, imgSpeedFout, imgZFout] = getAddr
            d = rdir('traceData*.mat');
        [Y,I] = sort([d.datenum]); 
        traceDataPath = {d(I(end)).name}; % load the most recent traceData
        %display(traceDataPath);
        % input file name
            imgFin = 'maxProj.tif'; % input : projection movie
        % output file name
        imgSpeedFout = 'Dorsal3D-wTracesCinSpeed.tif'; % output : color in Speed
        imgZFout = 'Dorsal3D-wTracesCinZ.tif'; % output : color in 
    end



    function plotTraces3D
        %% A - 3D data 

        %% 1- SET INPUT & OUTPUT
%        getPaths; % file paths
        %% 2- remove nan data and short traces
        FrstFrm = TraceX(:,1);
        TraceX = TraceX(~isnan(FrstFrm),:);
        TraceY = TraceY(~isnan(FrstFrm),:);
        TraceZ = TraceZ(~isnan(FrstFrm),:);
        zMax    = max(TraceZ(:));

        nzNumel = TraceX > 1;
        nzNumel = sum(nzNumel,2);
        nzNumel = nzNumel > minLenTrace;
        TraceX = TraceX(nzNumel,:);
        TraceY = TraceY(nzNumel,:);
        TraceZ = TraceZ(nzNumel,:);
        
        %% select one of the traces
        selTraces = 1;
        if selTraces
            TraceX2 = TraceX(R1:R2,:);
            TraceY2 = TraceY(R1:R2,:);
            TraceZ2 = TraceZ(R1:R2,:);
        end
        %% 3- CALC SPEED
        x = logical(TraceX2);
        x2 = logical(zeros(size(x)));
        x2(:,2:end) = x(:,1:end-1);
        x3= logical((x2).*x); % logical array refering to px where speed can be evaluated

        TraceXlast = zeros(size(TraceX2));
        TraceXlast(:,2:end) = TraceX2(:,1:end-1);
        TraceXdiff = TraceX2-TraceXlast;
        TraceXdiff = TraceXdiff.*(x3);

        TraceYlast = zeros(size(TraceY2));
        TraceYlast(:,2:end) = TraceY2(:,1:end-1);
        TraceYdiff = TraceY2-TraceYlast;
        TraceYdiff = TraceYdiff.*(x3);

        TraceZlast = zeros(size(TraceZ2));
        TraceZlast(:,2:end) = TraceZ2(:,1:end-1);
        TraceZdiff = TraceZ2-TraceZlast;
        TraceZdiff = TraceZdiff.*(x3);

        TraceDiff = sqrt(double(TraceXdiff.^2 + TraceYdiff.^2 + TraceZdiff.^2))*pxSz/tAcq;

        szNtr = size(TraceDiff,1);
        TraceDiffsmooth = zeros(size(TraceDiff));
        for i = 1:szNtr
            ixNZ=find(TraceDiff(i,:)~=0); % index of nonzero data
            TraceDiffsmooth(i,ixNZ) = smooth(TraceDiff(i,ixNZ),smthWin); % only smooth nonzero values
        end

        if smoothing
            TraceDiff = TraceDiffsmooth;
        end

        %% 4- SET figure properties
        imageInfo = imfinfo(imgFin);
        imSz=[imageInfo(1).Width,imageInfo(1).Height];
        gr = 0.5; % gray level
        GR = [gr gr gr];
        BL = [0 0 1];
        set(gca,'color',GR);
        colormap(CM);
        xlim([1 imSz(1)]);
        ylim([1 imSz(2)]);
        magImg = 2; % image magnification
        figSz(1) = imSz(1)*magImg;
        figSz(2) = imSz(2)*magImg;
            set(gca,'position',[0 0 1 1]); set(gcf,'position',[0 0 figSz(1) figSz(2)]);

        %% 5- CONFIG FIG DATA
        nonZeroSpeedIx = find(TraceDiff>0);
        
        if isSpeedinLog == 1
            TraceDiff = log(TraceDiff);
        end
        if isZinColor
            if isZinColor3D
                Cplot = TraceZ(nonZeroSpeedIx);
                Cplot = Cplot*pxSzZ;
                Zplot = TraceZ(nonZeroSpeedIx);   
            else
                Cplot = TraceZ(nonZeroSpeedIx);
                Cplot = Cplot*pxSzZ;
                Zplot = TraceDiff(nonZeroSpeedIx);         
            end
        else %speed in color
            Cplot = TraceDiff(nonZeroSpeedIx);
            %Cplot(Cplot>100)=0;
            %Cplot(Cplot<20)=0; % mu1-2
            if upLev 
                Cplot(find(Cplot>upLev)) = upLev;
            end
            if lowLev
                Cplot(find(Cplot<lowLev)) = lowLev;
            end
            Zplot = TraceZ(nonZeroSpeedIx);
        end
        maximum= max(max(Cplot));
        minimum= min(min(Cplot(find(Cplot>0))));
        
        %maximum = 20;
        %minimum = 50;
        %figure; hist(Cplot,20)
        figure(hFig);
        %% 6- PLOT
        plot3k({TraceX(nonZeroSpeedIx) imSz(2)-TraceY(nonZeroSpeedIx) Zplot},Cplot, [minimum maximum], dispPx);
        
        fprintf('colorbar: \n minimum= %d \n maximum = %d \n',minimum, maximum);

        %view(0,0); % XZ view
        axis equal;  
        axis([0 imSz(1) 0 imSz(2)]); 
        axis vis3d
        view(90,-90); 
        
        return;
        %% 7- EXPORT figure
        if isZinColor
            export_fig(imgZFout)
        else
            export_fig(imgSpeedFout)
        end    
    end

    %% B- plot version2
    function movieTraces3D
        hQ = 0;
        for ixFrm = 2:numFrames
            img2D = imread(imgFin,ixFrm); 
            delete(hImg);
            hImg = image(64*img2D/maxVal); axis image;
            set(gca,'position',[0 0 1 1]); set(gcf,'position',[0 0 figSz(1) figSz(2)]);
            ixTrc = find(TraceX(:,ixFrm)>0);
            currX = TraceX(:,ixFrm);
            currY = TraceY(:,ixFrm);
            currZ = TraceZ(:,ixFrm);
            uistack(hImg,'bottom');
            ixTrc2 = find(lastX(ixTrc)>0);
            ix = ixTrc(ixTrc2); % indices of traces tracked in the current frame    
            dspTrcIx = find(ixFrm == dspTrcFrm); % index for dissappearing traces
            hQdel = hQ(dspTrcIx,:); % handles for discont. traces
            delete(hQdel(find(hQdel~=0)));    % remove the traces of the discontinued traces
            for i = 1:length(ix) % for each trace 
                iL = ix(i);  % index for each line
                hQ(iL,ixFrm-1)=quiver(lastX(iL),lastY(iL),currX(iL)-lastX(iL),currY(iL)-lastY(iL),0,'Color',CM(round(currZ(iL)/zMax*64),:));
            end
            if ixFrm == 2
                export_fig(imgZFout)
            else
                export_fig(imgZFout,'-append');
            end
            fprintf('frame %i/%i \n',ixFrm,numFrames);
            lastX = currX; 
            lastY = currY;
            lastZ = currZ;
        end

        hold off; % release the figure
    end

end
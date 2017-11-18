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
    rotAng = 18; % Bi-Chang
    pxSzZ= 0.4; %[um]
    pxSzZ= pxSzZ/cos(rotAng*pi/180); %um
    %pxSzZ= 0.240; %um
    %% parameters
    isRunAll = 0;
    isRunWhole = 0; % not separated into D&V
    %% display param
    isDispAllXZ = 0;
    isDispAllXY = 1;
    isDispAllinOneFig = 0;
    isZinColor = 0; % ow speed in color
    isZinColor3D = 0;
    isTracesInColor = 0;
    dispPx = {'o',1}; % plot3k
    CM = ('jet'); %('jet');
    isSpeedinLog = 0;
    upLev = 0; % no speed color leveling
    upLev = 6*pxSz/tAcq;
    lowLev = 0; 
    lowLev = 3/8;
    vibrCorr = 0; % vibr. corrected
    smthWin = 5; % smoothing
    smoothing = 1;
    minLenTrace = 15;
    DVoverlap = 4;
    
    if isRunAll
        % FIND DATASET folders
        ixRt =1; goRoot; 
        CD = cd;
        %% get the folder names of DS(raw) and 3Drotated data
        if ismac 
            fPathDorsal = rdir('data/Rotated3D/**/*dorsal3D','isdir==1');
            fPathVentral = rdir('data/Rotated3D/**/*ventral3D','isdir==1');
        else % PC
            fPath3D = rdir('data\Rotated3D\**\*3D','strfind(name, ''\3D'')','isdir==1');
            fPathDorsal = rdir('data\Rotated3D\**\*dorsal3D','isdir==1');
            fPathVentral = rdir('data\Rotated3D\**\*ventral3D','isdir==1');
        end
        TraceX = [];
        TraceY = [];
        TraceZ = [];
        TraceINT = [];
        TraceT0 = [];
        nData = size(fPathDorsal,1); % number of datasets
        %fRange1 =1; fRange2=nData;
        fRange1 =   5; fRange2=5;
        hFig = figure; 
        nPlots = fRange2-fRange1+1;% number of plots
    else
        [traceDataPath, imgFin, imgSpeedFout, imgZFout] = getAddr;
        %traceDataPath = 'traceData.mat'; % vibration corrected
        
        load(char(traceDataPath)); % load traceData: 'TraceX','TraceY','TraceZ','TraceINT','TraceT0'
        % select traces
        hFig=figure;
        plotTraces3D
        saveFigImg('outImage.tif',999,hFig)
        return;
    end
    for i = fRange1 : fRange2 % runAll
        %% DS folder (rotAngle & phiTheta)
        goRoot;    
        %% 3D whole cell
        if isRunWhole
            cd(fPath3D(i).name);
            fprintf('Data folder: %s \n',fPath3D(i).name);
            [traceDataPath, imgFin, imgSpeedFout, imgZFout] = getAddr;
            load(char(traceDataPath)); % load traceData: 'TraceX','TraceY','TraceZ','TraceINT','TraceT0'
            % check # px in Z
            imageInfo = imfinfo('stack_001.tif');
            dorsalHeight = length([imageInfo.FileSize]); % # of layers
            if ceil(max(TraceZ(:))) ~= dorsalHeight
                fprintf('WARNING: dorsal height is not correct maxZ: %d dorsalHeight: %d \n', ceil(max(TraceZ(:))),dorsalHeight)
            end
        else
            %% Ventral
            cd(fPathVentral(i).name);
            fprintf('Data folder: %s \n',fPathVentral(i).name)
            [traceDataPath, imgFin, imgSpeedFout, imgZFout] = getAddr;
            load(char(traceDataPath)); % load traceData: 'TraceX','TraceY','TraceZ','TraceINT','TraceT0'
            if isDispAllinOneFig
                hFig=subplot(3,nPlots,(i-fRange1)+1)
                plotTraces3D
            else
                hFig=figure;
            end
            % check # px in Z
            imageInfo = imfinfo('stack_001.tif');
            ventralHeight = length([imageInfo.FileSize]); % # of layers
            TraceXv = TraceX; % ventral 
            TraceYv = TraceY;
            TraceZv = TraceZ;
            if ceil(max(TraceZ(:))) ~= ventralHeight
                fprintf('ERROR: ventral height is not correct maxZ: %d ventralHeight: %d \n', ceil(max(TraceZ(:))),ventralHeight)
            end
            goRoot;
            %% Dorsal
            cd(fPathDorsal(i).name);
            [traceDataPath, imgFin, imgSpeedFout, imgZFout] = getAddr;
            load(char(traceDataPath)); % load traceData: 'TraceX','TraceY','TraceZ','TraceINT','TraceT0'
            if isDispAllinOneFig
                subplot(3,nPlots,(i-fRange1)+1+nPlots)
                plotTraces3D
            end
            % check # px in Z
            imageInfo = imfinfo('stack_001.tif');
            dorsalHeight = length([imageInfo.FileSize]); % # of layers
            if ceil(max(TraceZ(:))) ~= dorsalHeight
                fprintf('WARNING: dorsal height is not correct maxZ: %d dorsalHeight: %d \n', ceil(max(TraceZ(:))),dorsalHeight)
            end
            %% fPathMain
            TraceX = vertcat(TraceXv, TraceX); % dorsal ventral combined
            TraceY = vertcat(TraceYv, TraceY);
            TraceZ = vertcat(TraceZv, TraceZ + ventralHeight-DVoverlap);
            clear TraceZv TraceYv TraceXv;
            if isDispAllinOneFig
                subplot(3,nPlots,(i-fRange1)+1+2*nPlots)
            end
        end
        %figure(i)
        plotTraces3D
        cd ..
        %return;
        view(-180,0) %CK : slide rotated z-view
        
   set(gca,'PlotBoxAspectRatioMode','auto'); set(gca,'CameraViewAngleMode','auto');
        return
        if ~isDispAllinOneFig
            if isZinColor
                %saveFigImg(strcat(sprintf('zImage_colorLev%0.2f-%0.2f',lowLevConst,upLevConst),'.tif'),999,hFig)
                saveFigImg(strcat(sprintf('zImage'),'.tif'),999,hFig)
            else
                saveFigImg(strcat(sprintf('speedImg_colorLev%0.2f-%0.2f',lowLev,upLev),'.tif'),999,hFig)
            end
        end
        
%         %% figure annotation 
%         if isZinColor 
%             titleTx = 'Z in color, speed in 3D ';
%         else
%             if isSpeedinLog
%                 titleTx = 'Speed (log) in color (3D)';
%             else
%                 titleTx = 'Speed in color (3D)';
%             end
%         end
%         ax = axes('position',[0,0,1,1],'visible','off');
%         tx = text(0.5,0.95,titleTx);
%         set(tx,'fontweight','bold');
%         tx = text(0.05,0.8,'Ventral'); 
%         set(tx,'fontweight','bold');
%         tx = text(0.05,0.5,'Dorsal'); 
%         set(tx,'fontweight','bold');
%         tx = text(0.05,0.2,'3D');
%         set(tx,'fontweight','bold');
        
        
    end


    function [traceDataPath, imgFin, imgSpeedFout, imgZFout] = getAddr

        [Y,I] = sort([d.datenum]); 
        traceDataPath = {d(I(end)).name}; % load the most recent traceData
        %display(traceDataPath);
        % input file name
        if strfind(cd,'ventral')
            imgFin = 'maxProjVentral3D.tif'; % input : projection movie
        elseif strfind(cd,'dorsal')
            imgFin = 'maxProjDorsal3D.tif'; % input : projection movie
        else
            imgFinAddr = dir('maxProj*')
            if isempty(imgFinAddr)
                disp('ERROR : no projection image is found')
            else
                imgFin = imgFinAddr(1).name;
            end
            
        end
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
        selTraces = 0;
        if selTraces
            TraceX = TraceX(5,:);
            TraceY = TraceY(5,:);
            TraceZ = TraceZ(5,:);
        end
        %% 3- CALC SPEED
        x = logical(TraceX);
        x2 = logical(zeros(size(x)));
        x2(:,2:end) = x(:,1:end-1);
        x3= logical((x2).*x); % logical array refering to px where speed can be evaluated

        TraceXlast = zeros(size(TraceX));
        TraceXlast(:,2:end) = TraceX(:,1:end-1);
        TraceXdiff = TraceX-TraceXlast;
        TraceXdiff = TraceXdiff.*(x3);

        TraceYlast = zeros(size(TraceY));
        TraceYlast(:,2:end) = TraceY(:,1:end-1);
        TraceYdiff = TraceY-TraceYlast;
        TraceYdiff = TraceYdiff.*(x3);

        TraceZlast = zeros(size(TraceZ));
        TraceZlast(:,2:end) = TraceZ(:,1:end-1);
        TraceZdiff = TraceZ-TraceZlast;
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
        if isDispAllXY || isDispAllXZ
            %set(gcf,'position',[0 0 figSz(1) figSz(2)]);
            %set(hFig,'WindowStyle','docked');
        else
            set(gca,'position',[0 0 1 1]); set(gcf,'position',[0 0 figSz(1) figSz(2)]);
        end

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
        elseif isTracesInColor
            CplotVecN = size(TraceX,1); % # of traces
            Nframe = size(TraceX,1); % # of frames
            CplotVec = 1: CplotVecN;
            Cplot = repmat(CplotVec',1,Nframe);
            Cplot = Cplot(nonZeroSpeedIx);
            Zplot = ones(size(TraceZ));  
            Zplot = Zplot(nonZeroSpeedIx);
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

        if isDispAllXZ
            view(0,0); % XZ view
        elseif isDispAllXY
            axis equal;  
            axis([0 imSz(1) 0 imSz(2)]); 
            axis vis3d
            view(90,-90); 
        end
        
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
clear all; close all;
    
    %load('traceData-coeff72')
    load('traceData2-coeff30')
    %load('strainData')
    disp('loading plot data!!!')
    % remove nan rows
    
        TraceX = full(TraceX);
        TraceY = full(TraceY);
        TraceZ = full(TraceZ);
        TraceSpeed = full(TraceSpeed);
        
        
    if 1 % thinsulate
        ix = 1:10:size(TraceX,1);
        TraceX = TraceX(ix,:);
        TraceY = TraceY(ix,:);
        TraceZ = TraceZ(ix,:);
        TraceSpeed = TraceSpeed(ix,:);
    end
        
    if 0
        TraceX = TraceX(find(~isnan(min(TraceX,[],2))),:);
        TraceY = TraceY(find(~isnan(min(TraceX,[],2))),:);
        TraceZ = TraceZ(find(~isnan(min(TraceX,[],2))),:);    
        TraceSpeed = TraceSpeed(find(~isnan(min(TraceX,[],2))),:);    
    end
    
    inputInfo = dir('inputInfo*.mat');
    load(inputInfo.name);
    pxSzX = 2.501;
    pxSzY = 2.501;
    pxSzZ = 6;
    pxSzX = 1;
    pxSzY = 1;
    pxSzZ = 1;
    pxszZvsXY  = (pxSzZ/pxSzX); % ratio of pixel sizes in XY and Z
    
    img = imread('maxProjDorsal3D.tif',1);
    [Boy1, En1] = size(img);
    
    
    if 0
        % clear traces on the edge
        TraceX = full(TraceX);
        TraceY = full(TraceY);
        TraceX(TraceX==0) = nan;
        TraceY(TraceY==0) = nan;

        mxX = max(TraceX,[],2);
        mnX = min(TraceX,[],2);
        mxY = max(TraceY,[],2);
        mnY = min(TraceY,[],2);
        discardTraces = find((mnX<=5) + (mxX>En1-4) + (mnY<=5) + (mnY>Boy1-4));
        for i = 1:numel(discardTraces)
            TraceX=TraceX([1:discardTraces(i)-1 discardTraces(i)+1:end],:);
            TraceY=TraceY([1:discardTraces(i)-1 discardTraces(i)+1:end],:);
            TraceZ=TraceZ([1:discardTraces(i)-1 discardTraces(i)+1:end],:);
            TraceSpeed=TraceSpeed([1:discardTraces(i)-1 discardTraces(i)+1:end],:);
            discardTraces = discardTraces - 1;
        end    
    end
    
    %Boy1 = 1024;
    isCropData = 0; 
    if isCropData
        Tr1=1;
        Tr2=500;
        %Tr2=size(TraceY,1);
        tt1 = 1; 
        tt2 = size(TraceY,2);
        tt2 = 50;
        dTr = 100;
        TraceX = TraceX(Tr1:dTr:Tr2,tt1:tt2);
        TraceY = TraceY(Tr1:dTr:Tr2,tt1:tt2);
        TraceZ = TraceZ(Tr1:dTr:Tr2,tt1:tt2);
        Strain = Strain(Tr1:dTr:Tr2,tt1:tt2);
        TraceSpeed = TraceSpeed(Tr1:dTr:Tr2,tt1:tt2);
    end
    
    if 0
        for i = 1:size(TraceX,1)
            X = full(mean(TraceX(i,TraceX(i,:)>0),2));
            Y = full(mean(TraceY(i,TraceX(i,:)>0),2));
            disp(sprintf('X:%.02f,Y:%.02f\n',X,Y));
        end
    end
        
    % color data
    CplotVecN = size(TraceX,1); % # of traces
    Nframe = size(TraceX,2); % # of frames
    CMinZ = 1;
    Sdata = TraceZ;
    %Sdata = Strain;
    isSpeed = 0; %Plots speed io strain
    if isSpeed, Sdata = TraceSpeed; end;
    
    if CMinZ
        CplotVec = 1:CplotVecN;
        CplotVec = rem(CplotVec,64)+1;
        Cplot = repmat(CplotVec',1,Nframe);
        CM = colormap('lines');
    else
        Cdata = Sdata;
        Cmin = min(Cdata(Cdata~=0));
        Cmax = max(Cdata(:));
        Crange = Cmax - Cmin;
        Cplot = round(63*(Cdata-Cmin) / Crange) + 1;
        CM = colormap('jet');
    end
    
    
    Frames = size(TraceX,2);
    isCropTr=0;
    if isCropTr
        Ntr0 = 1; Ntr = size(TraceX,1);
        trVec = [Ntr0:Ntr];
        TraceX = TraceX(trVec,1:Frames);
        TraceY = TraceY(trVec,1:Frames);
        TraceZ = TraceZ(trVec,1:Frames);
    end

    En = Frames;
    [Boy2]=size(TraceX,1);
    
    nonZeroIx = find(TraceY>0);
    dispPx = {'*',2}; % plot3k
    minimum = 0;
    maximum = max(Cplot(:));
    fig1=figure;
    colormap('lines')

    hold on

    subplot(2,5,[1:4 6:9]);
        xx = TraceX(TraceSpeed>0);
        yy = TraceY(TraceSpeed>0);
        zz = TraceZ(TraceSpeed>0);
        cc = CM(Cplot(TraceSpeed>0),:);
    scatter3(xx,Boy1-yy,zz,4,cc,'*');
    %plot3k({TraceX(nonZeroIx) Boy1-TraceY(nonZeroIx) pxszZvsXY*TraceZ(nonZeroIx)},Cplot(nonZeroIx), [minimum maximum], dispPx);
    view(45,20);
    view(0,90);
    xlim([0 En1]); ylim([0 Boy1]);
    hold off
    htitle = title('traces in 3D')
    axis equal
    
    %% plot selected traces
    TraceX(TraceX==0) = nan;
    TraceY(TraceX==0) = nan;
    maximize
    % frame number input
    %hEdit = uicontrol('style','edit','BackgroundColor',[1 1 1]);
    
    while 1
        isSel = 0;
        while ~isSel
            [x,y] = ginput(1);
            y = Boy1-y;
            disp(sprintf('X:%.02f,Y:%.02f\n',x,y));
            dist = sqrt((TraceX-x).^2 + (TraceY-y).^2);
            [val ix] = min(dist(:));
            if val < 5
                isSel = 1;
            end
            if x>-200 && x<-63 && y>975 && y< 992
                isSel = 1;
                isStop = 1;
            end
            
        end
        
        % trace number
        [tr fr] = ind2sub(size(TraceX),ix);
        traceNumInput = get(hEdit,'String');
        if ~isempty(traceNumInput), tr = str2num(traceNumInput); end;
            
        set(htitle,'String',sprintf('trace#:%i',tr))
        fr = find(TraceX(tr,:)>0);
        subplot(2,5,[5]);
        tx = TraceX(tr,fr); mtx = min(tx);
        tx = tx - mtx;
        ty = TraceY(tr,fr); mty = min(ty);
        ty = ty - mty;
        tz = TraceZ(tr,fr); mtz = min(tz);
        tz = tz - mtz;
        ts = TraceSpeed(tr,fr);
        sDataText = 'speed';
        if isSpeed, sDataText = 'TraceSpeed'; end;
        [hAx,hLine1,hLine2] = plotyy(fr,tx,fr,ts); hold on;
        set(hLine1,'Color',[0 0 1])
        set(hLine2,'Color',[0 0 0])
        hAx(3) = plot(fr,ty,'g-');
        hAx(4) = plot(fr,tz,'r-');
        legend('X','Y','Z',sDataText);
        ylabel(hAx(1),'position [px]') % left y-axis
        ylabel(hAx(2),[sDataText ' speed [px/frame]']) % right y-axis
        xlabel('frames')
        for i = 1:2, xlim(hAx(i),[0 Frames+2] ); end;
        v = [tx ty tz];
        ylim(hAx(1),[min(v)-2 max(v)+2] )
        ylim(hAx(2),[min(ts) max(ts)] )
%        ylim(hAx(3),[min(v)-2 max(v)+2] )
%        ylim(hAx(4),[min(v)-2 max(v)+2] )
        title([sDataText '- XYZ'])
        hold off;

        subplot(2,5,[10]);
        plot(TraceX(tr,fr),Boy1-TraceY(tr,fr));
        xlabel('X')
        ylabel('Y') % right y-axis
        title('XY trajectory')
        axis equal

    end

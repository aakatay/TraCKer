clear all; close all;
    load('traceData-coeff96')
    load('strainData')
    disp('loading plot data!!!')
    
    % remove nan rows
    TraceX = TraceX(find(~isnan(min(TraceX,[],2))),:);
    TraceY = TraceY(find(~isnan(min(TraceX,[],2))),:);
    TraceZ = TraceZ(find(~isnan(min(TraceX,[],2))),:);    
    
    pxSzX = 2.501;
    pxSzY = 2.501;
    pxSzZ = 6;
    pxszZvsXY  = (pxSzZ/pxSzX); % ratio of pixel sizes in XY and Z
    
    Boy1 = 1024;
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
        Speed = Speed(Tr1:dTr:Tr2,tt1:tt2);
    end
    
    for i = 1:size(TraceX,1)
        X = mean(TraceX(i,TraceX(i,:)>0),2);
        Y = mean(TraceY(i,TraceX(i,:)>0),2);
        disp(sprintf('X:%.02f,Y:%.02f\n',X,Y));
    end
        
    % color data
    CplotVecN = size(TraceX,1); % # of traces
    Nframe = size(TraceX,2); % # of frames
    CMinZ = 0;
    if CMinZ
        CplotVec = 1: CplotVecN;
        CplotVec = rem(CplotVec,64)+1;
        Cplot = repmat(CplotVec',1,Nframe);
    else
        Cdata = Strain;
        Cmin = min(Cdata(Cdata~=0));
        Cmax = max(Cdata(:));
        Crange = Cmax - Cmin;
        Cplot = round(63*(Cdata-Cmin) / Crange) + 1;
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
    minimum = 1;
    maximum = max(Cplot(:));
    fig1=figure;
    colormap('lines')

    hold on

    subplot(2,5,[1:4 6:9]);
    plot3k({TraceX(nonZeroIx) Boy1-TraceY(nonZeroIx) pxszZvsXY*TraceZ(nonZeroIx)},Cplot(nonZeroIx), [minimum maximum], dispPx);
    view(45,20);
    view(0,90);
    xlim([0 Boy1]); ylim([0 Boy1]);
    hold off
    title('traces in 3D')
    
    %% plot selected traces
    TraceX(TraceX==0) = nan;
    TraceY(TraceX==0) = nan;
    maximize
    
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

        [tr fr] = ind2sub(size(TraceX),ix);
        fr = find(TraceX(tr,:)>0);
        subplot(2,5,[5]);
        tx = TraceX(tr,fr); mtx = min(tx);
        tx = tx - mtx;
        ty = TraceY(tr,fr); mty = min(ty);
        ty = ty - mty;
        tz = TraceZ(tr,fr); mtz = min(tz);
        tz = tz - mtz;
        ts = Strain(tr,fr);
        [hAx,hLine1,hLine2] = plotyy(fr,tx,fr,ts); hold on;
        set(hLine1,'Color',[0 0 1])
        set(hLine2,'Color',[0 0 0])
        plot(fr,ty,'g-');
        plot(fr,tz,'r-');
        legend('X','Y','Z','strain');
        ylabel(hAx(1),'position [px]') % left y-axis
        ylabel(hAx(2),'strain [px/frame]') % right y-axis
        xlabel('frames')
        %xlim([1 Frames] )
        title('strain - XYZ')
        hold off;

        subplot(2,5,[10]);
        plot(fr,Speed(tr,fr));
        xlabel('frames')
        ylabel('speed [px/frame]') % right y-axis
        title('speed')

    end

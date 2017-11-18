
    close all;
    F = findall(0,'type','figure','tag','TMWWaitbar'); delete(F);
    pxSzX = 2.501;
    pxSzY = 2.501;
    pxSzZ = 6;
    sptReAppearTime = 2; % use the same value used in the TraCKer
    magImg = 6; % plots magnified image
    pxszZvsXY  = (pxSzZ/pxSzX); % ratio of pixel sizes in XY and Z
    df = 3; % # the number of frames between two frames used for strain calc.
    
    evalWinRpx = 85; % [px] evaluation window radius
    isSelFirst5Neighbour = 0; % uses the closest 5 neighbour traces
    evalWinR = evalWinRpx*pxSzX; 

    d = rdir('traceData-coeff*.mat');
    [Y,J] = sort([d.datenum]); 
    traceDataPath = {d(J(end)).name}; % load the most recent posData-coeff
    load(char(traceDataPath));
    X = TraceX(find(~isnan(min(TraceX,[],2))),:);
    Y = TraceY(find(~isnan(min(TraceX,[],2))),:);
    Z = TraceZ(find(~isnan(min(TraceX,[],2))),:);    
    INT = TraceINT(find(~isnan(TraceX(:,1))),:);
    clear TraceX TraceY TraceZ

    [traces frames] = size(X);
    %% crop
    isCrop = 0;
    if isCrop
        t1 = 174; t2 = 177;
        X = X(t1:t2,:);
        Y = Y(t1:t2,:);
    else
        t1 = 1;
        t2 = traces;
    end
    %% calculate displacement
    for tr = 1:traces
        %% fill gaps (missing frames in jumpy traces)
        traceX = X(tr,:);
        traceY = Y(tr,:);
        traceZ = Z(tr,:);
        ind = find(traceX~=0);
        frst = min(ind); last = max(ind);
        %iJump = find(trace(frst:last)==0)+frst-1;

        jmp = abs((traceX(frst:last)>0)-1);
        if ~sum(jmp > 0)
            TraceX(tr,:) = X(tr,:);
            TraceY(tr,:) = Y(tr,:);
            TraceZ(tr,:) = Z(tr,:);
        else
            bnd = bwboundaries(jmp,'noholes');
            for i = 1:numel(bnd) % for each jump
                temp = bnd{i}+frst-1; % jump boundaries
                jb = temp(:,2);
                mx= max(jb); mn=min(jb);
                jL = mx-mn+2; % length

                jSx = traceX(mx+1)-traceX(mn-1);% size
                jSy = traceY(mx+1)-traceY(mn-1);% size
                jSz = traceZ(mx+1)-traceZ(mn-1);% size
                jsX = jSx/jL;% step
                jsY = jSy/jL;% step
                jsZ = jSz/jL;% step
                jVx = traceX(mn-1)+jsX*(1:jL-1);
                jVy = traceY(mn-1)+jsY*(1:jL-1);
                jVz = traceZ(mn-1)+jsZ*(1:jL-1);
                traceX(mn:mx)=jVx;
                traceY(mn:mx)=jVy;
                traceZ(mn:mx)=jVz;
                TraceX(tr,:) = traceX;
                TraceY(tr,:) = traceY;
                TraceZ(tr,:) = traceZ;
            end
        end
        traceXdiff = [0 traceX] - [traceX 0];
        traceXdiff = traceXdiff(2:end-1);
        traceYdiff = [0 traceY] - [traceY 0];
        traceYdiff = traceYdiff(2:end-1);
        traceZdiff = [0 traceZ] - [traceZ 0];
        traceZdiff = traceZdiff(2:end-1);
        traceDisp(tr,:) = sqrt(traceXdiff.^2 + traceYdiff.^2 + pxszZvsXY.^2 *traceZdiff.^2);
        if frst ~= 1
            traceDisp(tr,frst-1) = 0;
        end
        if last ~= frames
            traceDisp(tr,last) = 0;
        end
    end
    %% calculate strain
    Strain = nan(traces,frames);
    Speed = nan(traces,frames);
    numNeighbour = 0; nSpots = 0;
    for f = 1:frames-1
        
        for tr = t1:t2 % traces
            if ~TraceX(tr,f) || numel(traceDisp(find(traceDisp(:,f)>0),f)) ==1
                continue;
            end
            Speed(tr,f) = traceDisp(tr,f);
            if f >frames-df, continue; end;
            Xc = TraceX(tr,f); 
            Yc = TraceY(tr,f);
            Zc = TraceZ(tr,f);
            findTracesAtProx; % [TrProxIx dR]

            % only first 5
            if isSelFirst5Neighbour
                [v ix] =sort(dR);
                TrProxIx = TrProxIx(ix(1:5)); 
            end

            numNeighbour = numel(TrProxIx)+numNeighbour ;
            nSpots = nSpots + 1;

            diffX = TraceX(TrProxIx2,f) - Xc;
            diffY = TraceY(TrProxIx2,f) - Yc;
            diffZ = TraceZ(TrProxIx2,f) - Zc;
            diff = sqrt(diffX.^2 + diffY.^2 + pxszZvsXY.^2*diffZ.^2); 
            diffX2 = TraceX(TrProxIx2,f+df) - Xc;
            diffY2 = TraceY(TrProxIx2,f+df) - Yc;
            diffZ2 = TraceZ(TrProxIx2,f+df) - Zc;
            diff2 = sqrt(diffX2.^2 + diffY2.^2 + pxszZvsXY.^2 *diffZ2.^2);
            DIFF = diff2 - diff;
            nNeighbour(tr,f) = numel(TrProxIx2);
            DIFF = sum(DIFF)/nNeighbour(tr,f);

            Strain(tr,f) = DIFF;
            if Strain(tr,f) == 0
                disp('error');
            end
        end
    end
    save('strainData.mat','Strain','Speed')
    disp(sprintf('evalWinRpx: %i average neighbours: %.02f\n',evalWinRpx,numNeighbour/nSpots)); %return
    
    %% generate colormap
    nTicks = 6;
    figure;  colorbarData = plotColorBar(Strain/df,nTicks);
    print(sprintf('colorbar_R%i.bmp',evalWinRpx),'-dbitmap');
    %save('colorbarData',colorbarData);
    fc = fopen('colorData.txt','w');
    for i = 1:numel(colorbarData)
        fprintf(fc,'%.02f\n',colorbarData(i));
    end
    fclose(fc);
    
    %% image file
    if exist('fname.mat')
        load fname
        imgFin = fname;
    else
        imgFin = 'maxProj.tif';
    end
    if ~exist(imgFin), imgFin2 = rdir(sprintf('..\\%s',imgFin)); imgFin=imgFin2.name; end;
    I = imread(imgFin);    
    if exist('cropCoor.mat')
        load cropCoor; 
        [yNpx,xNpx]=size(I(yy1:yy2,xx1:xx2));
    else
        imageInfo = imfinfo(fname);
        imSz=[imageInfo(1).Width,imageInfo(1).Height];
        yy1=1;xx1=1; yy2 = imSz(2); xx2 = imSz(1);
        [yNpx,xNpx]=size(I);
    end;
    
    
    
%% pixelate the coordinates
    isPx = 0;
    if isPx
        szGaus = 100; sgGaus = 70;
        gaus=fspecial('gaussian', szGaus, sgGaus);
        figure; maximize
        for t = 1:frames % # of frames
            % pixelate the coordinates
            x = TraceX(~isnan(Strain(:,t)),t);
            y = TraceY(~isnan(Strain(:,t)),t);
            x = ceil(x);
            y = ceil(y);
            A = single(zeros(yNpx,xNpx));
            A(sub2ind([yNpx,xNpx],y,x)) = Strain(~isnan(Strain(:,t)),t); % pits
            A2(:,:,t) = A;
            imagesc(A)
            waitforbuttonpress
        end
        %close(hh)
    end
    
    %% contour plot
    demag = 8;
    szContImg = yNpx/demag;
    if ~exist('pxStrain.mat')
        contImg = zeros(szContImg,szContImg);
        xx = 1:szContImg;
        [XX YY] = meshgrid(xx,xx);
        XX1 = XX*demag; YY1 = YY*demag;
        xx = 1:yNpx;
        [XX2 YY2] = meshgrid(xx,xx);
        hh = waitbar(0,'generating strain map...');
        pxStrain = zeros(yNpx,yNpx,frames-df);
        for f = 1:frames-df
            for i = 1:szContImg
                for j = 1:szContImg
                    Xc = i*demag-demag/2; Yc = j*demag-demag/2;
                    findTracesAtProx; % TrProxIx2
                    pxStrain0(szContImg-j+1,i) = sum(Strain(TrProxIx2,f))/numel(TrProxIx2); % pixelated strain
                end
            end 
            temp = interp2(XX1,YY1,pxStrain0,XX2,YY2);
            pxStrain(demag/2+1:yNpx-demag/2,demag/2+1:yNpx-demag/2,f) = temp(demag+1:yNpx,demag+1:yNpx);
            waitbar(f/(frames-df))
        end
        save('pxStrain','pxStrain')
        close(hh)
    else
        load('pxStrain')
    end
        
        
    %% PLOT
    % color data
    CplotVecN = size(TraceX,1); % # of traces
    Nframe = size(TraceX,2); % # of frames
    CMinZ = 0;
    Cdisp = 1; % coloring by displacement
    if CMinZ
        CplotVec = 1: CplotVecN;
        CplotVec = rem(CplotVec,64)+1;
        Cplot = repmat(CplotVec',1,Nframe);
    else
        if Cdisp         
            Cdata = Strain; % strain
        else
           Cdata =  TraceZ;
        end
        Cmin = min(Cdata(Cdata~=0));
        Cmax = max(Cdata(:));
        Crange = Cmax - Cmin;
        Cplot = round(63*(Cdata-Cmin) / Crange) + 1;
    end

    Frames = frames;
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
    
%%
    NAN = find(isnan(TraceX));
    TraceX(NAN)=0;
    TraceY(NAN)=0;
    TraceZ(NAN)=0;

    % parameters:
    %load(posDataFileNm,'sptReAppearTime'); %(frames) use the value from tracker function generating trace values
    isTrace = 0;
    isStrain = 1;

    load fname;
    imgZFout = 'TraceImage';
    TraceY = yNpx-TraceY+1;

    figSz(1) = xNpx*magImg;
    figSz(2) = yNpx*magImg;
    CM = colormap('jet');
  

    % find the frames where the traces disappear
    hat = ones(1,sptReAppearTime); 
    dspTrcFrm = zeros(size(TraceX,1),1);
    for i = 1:size(TraceX,1) % all traces
        trFrm = conv(TraceX(i,:),hat,'same'); % active frames of the trace
        trFrm(1:find(trFrm>0,1)) = 1;
        lastFrm = find(trFrm==0,1); 
        if ~isempty(lastFrm)
            dspTrcFrm(i) = lastFrm;
        end
    end

    hQ = 0; %hImg = image; 
    lastX = TraceX(:,1); 
    lastY = TraceY(:,1);
    frm1 = 1;
    tit = 'image';
    m = xNpx;
    n = yNpx;
    pos=get(0,'ScreenSize');
    pos=pos(3:4) - [m n-35];
    
    isSurf = 1;
    if isSurf, colormap('jet'); else, colormap('gray'); end;
    for ixFrm = frm1:En-df
        fig=figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Position',[pos/2 m n]);
        axe=axes('Parent',fig,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren','XLim',0.5+[0 m],'YLim',0.5+[0 n]);
    
        img2D = imread(fname,ixFrm); 
        img2D = flipud(img2D(yy1:yy2,xx1:xx2));
        %delete(hImg);
        if isSurf
            hCont = contour(pxStrain(demag/2+1:yNpx-demag/2,demag/2+1:yNpx-demag/2,ixFrm));
            xlim([-demag/2+1 yNpx-demag/2]); ylim([-demag/2+1 yNpx-demag/2]); 
        else
            hImg = imagesc(img2D,'Parent',axe); %axis image; 
        end
        %set(gca,'position',[0 0 1 1]); set(gcf,'position',[0 0 figSz(1) figSz(2)]);
        ixTrc = find(TraceX(:,ixFrm)>0);
        currX = TraceX(:,ixFrm);
        currY = TraceY(:,ixFrm);
        hold on
        ixTrc2 = find(lastX(ixTrc)>0);
        if isStrain
            ix = find(Cplot(:,ixFrm)>0);
        else
            ix = ixTrc(ixTrc2); % indices of traces tracked in the current frame    
        end
        dspTrcIx = find(ixFrm == dspTrcFrm); % index for dissappearing traces
        showTrace =1;
        if showTrace && isTrace
            hQdel = hQ(dspTrcIx,:); % handles for discont. traces
            delete(hQdel(find(hQdel~=0)));    % remove the traces of the discontinued traces
        end
        if isStrain            
            hScat = scatter(currX(ix),currY(ix),20,CM(Cplot(ix,ixFrm),:),'filled','o');
        elseif isTrace 
            for i = 1:round(length(ix)) % for each trace 
                iL = ix(i);  % index for each line
                hQ(iL,ixFrm-1)=quiver(lastX(iL),lastY(iL),currX(iL)-lastX(iL),currY(iL)-lastY(iL),'Color',CM(Cplot(iL,ixFrm),:));
                if ~showTrace
                    adjust_quiver_arrowhead_size(hQ(iL,ixFrm-1),5)
                end
            end
        end
        dig = floor(log10(Frames))+1;
        resLow = 120; % [dpi] resolution
        rLow = sprintf('-r%i',resLow);
        switch dig
            case 1
                print(sprintf('%s_%01i.tif',imgZFout,ixFrm),'-dtiff',rLow); 
            case 2
                print(sprintf('%s_%02i.tif',imgZFout,ixFrm),'-dtiff',rLow); 
            case 3
                print(sprintf('%s_%03i.tif',imgZFout,ixFrm),'-dtiff',rLow); 
            case 4
                print(sprintf('%s_%04i.tif',imgZFout,ixFrm),'-dtiff',rLow); 
        end
        if exist('hScat'), delete(hScat); end;
        %if isSurf, ; else, delete(hImg); end;
        hold off

        %waitforbuttonpress
        if ~showTrace
            delete(hQ(find(hQ(:,ixFrm-1)~=0),ixFrm-1));
        end

        %export_fig(imgZFout,'-append');
        fprintf('frame %i/%i \n',ixFrm,En);
        lastX = currX; 
        lastY = currY;
        close(fig)
%        lastZ = currZ;

    end
    tiff2stack('TraceImage'); % combine frames in a stack
    


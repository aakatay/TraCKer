
    clear;
    close all;
    F = findall(0,'type','figure','tag','TMWWaitbar'); delete(F)

    imgFin = 'maxProj.tif';
    if ~exist(imgFin), imgFin2 = rdir(sprintf('..\\%s',imgFin)); imgFin=imgFin2.name; end;
    I = imread(imgFin);
    [yNpx,xNpx]=size(I);
    imgMax= max(I(:));
    
    dT = 1; % time interval [sec]
    %% load xyz data
    isTraceData = 1;
    if isTraceData
        d = rdir('traceData0-coeff*.mat');
        [Y,J] = sort([d.datenum]); 
        traceDataPath = {d(J(end)).name}; % load the most recent posData-coeff
        load(char(traceDataPath));

        TraceX(isnan(TraceX)) = 0;
        TraceY(isnan(TraceY)) = 0;
        TraceZ(isnan(TraceZ)) = 0;
        TraceINT(isnan(TraceINT)) = 0;

        isGenData = 0;
        if isGenData
            genCellContractPits
        end            
        X = TraceX;
        Y = TraceY;
        Z = TraceZ;
        INT = TraceINT;
        INT(find(INT<0))=min(INT(find(INT>0)));
    else % posData
        d = rdir('xyzData-coeff*.mat');
        [Y,J] = sort([d.datenum]); 
        xyzDataPath = {d(J(end)).name}; % load the most recent posData-coeff
        load(char(xyzDataPath));
        INT(find(INT<0))=min(INT(find(INT>0)));
        %load('xyzData-coeff1979');
    end
    %% colormap
    colormap('jet'); 
    CM = colormap;
    intMax = max(INT(:));
    f = 1; % frame #
    R1 = 1; R2 = 10;
    C = CM(round(INT(R1:R2,f)/intMax/2*64+32),:);
    %scatter3(X(R1:R2,f),Y(R1:R2,f),Z(R1:R2,f),INT(R1:R2,f)/1000,C)
    
    %% # of spots at each frame
    T = size(X,2); % # of frames
    for t = 1:T
        nSpotsTime(t) = 1+0*numel(find(X(:,t)>0));
    end
    %plot(nSpotsTime);
    
    %% mask the data for one cell
    isMask = 0;
    if isMask
        load('outlineMask');
        X2 = zeros(size(X));
        Y2 = zeros(size(Y));
        for t = 1:T % # of frames
            % pixelate the coordinates
            x = X(find(X(:,t)~=0),t);
            y = Y(find(X(:,t)~=0),t);
            %z = Z(find(X~=0));
            x = ceil(x);
            y = ceil(y);
            inCell = BW(sub2ind(size(BW),x,y,t*ones(size(x))))==1;
            x = x.*inCell;
            y = y.*inCell;
            X2(1:size(x,1),t) = x;
            Y2(1:size(x,1),t) = y;
        end
        X = X2;
        Y = Y2;
    end

    fnameOut_ro = '_roStack.tif';
    fnameOut_pits = '_pitsStack.tif';
    szGaus = 100; sgGaus = 70;
    gaus=fspecial('gaussian', szGaus, sgGaus);
    figure; maximize
    h = waitbar(0,'calculating the density map...');
    for t = 1:T % # of frames  
        x = X(find(X(:,t)~=0),t);
        y = Y(find(X(:,t)~=0),t);
        I = INT(find(X(:,t)~=0),t);
        %z = Z(find(X~=0));
        x = ceil(x);
        y = ceil(y);
        A = single(zeros(yNpx,xNpx));
        A(sub2ind([yNpx,xNpx],y,x)) = 1; % pits
        neigbR = 25; % the neighbouring distance (radius [px])
        a = [-neigbR:neigbR];
        cX = meshgrid(a);
        cY = meshgrid(a)';
        cZ = neigbR^2-cX.^2-cY.^2;
        cZ(find(cZ<0))=0;
        convF = gaus-gaus(szGaus/2,1);
        convF(find(convF<0)) = 0;
        %A = zeros(size(A));
        %A(250,250) = 1; A(250,288) = 1;
        Aconv = filter2(convF,A);
        %Aconv = conv2(A,convF,'same');
        %imagesc(Aconv)
        roStack(:,:,t) = Aconv/nSpotsTime(t); % x,y,t stack keeping density map

        if t>1
            roStackdiff(:,:,t-1) = abs(roStack(:,:,t) - roStack(:,:,t-1));
            %imagesc(roStackdiff(:,:,t-1))
        end
        subplot(1,2,1)
        imagesc(roStack(:,:,t)); axis image;
        %hold;

        C = CM(ceil(I/intMax*64),:);
        %scatter(x,y,I/1000,C,'.')
        hold off;
        colormap('jet')
        colormap('autumn')
        colormap('copper')
        cm = colormap; cm(1,:)= [0 0 0];
        colormap(cm)
        pitStack(:,:,t) = Aconv.*A/nSpotsTime(t); % density of pits coded to pit positions
        subplot(1,2,2)
        imagesc(pitStack(:,:,t)); axis image;
        %export_fig('imgOut.tif');
        title(sprintf('time:%i',t*dT))
       %pause(0.15)
        waitbar(t/T);
    end

    savedata = roStack;
    savedata = 2^8*savedata/max(savedata(:));
    for t = 1:size(savedata,3) % # of frames  
        if t == 1
            imwrite(uint16(savedata(:,:,t)),fnameOut_ro);
        else
            imwrite(uint16(savedata(:,:,t)),fnameOut_ro,'tiff','WriteMode','append');
        end
    end
    savedata = pitStack;
    savedata = single(imgMax)*savedata/max(savedata(:));
    for t = 1:size(savedata,3) % # of frames  
        if t == 1
            imwrite(uint16(savedata(:,:,t)),fnameOut_pits);
        else
            imwrite(uint16(savedata(:,:,t)),fnameOut_pits,'tiff','WriteMode','append');
        end
    end    
    close(h)
    waitforbuttonpress
    close all;
    return

    [TrIx_,FrIx_]=ind2sub(size(TraceX),find(TraceX>0)); % trace and frame indices

for i = 1:numel(TrIx_)
    TrIx = TrIx_(i); FrIx = FrIx_(i);
    TraceRo(TrIx,FrIx) = roStack(ceil(TraceX(TrIx,FrIx)),ceil(TraceY(TrIx,FrIx)),FrIx); % density
end

for i = 1:nTrace
    meanTr(i) = mean(TraceRo(i,find(TraceRo(i,:))>0));
    TraceRo(i,:) = TraceRo(i,:)/meanTr(i);
end
pitStack2 = [];
for t = 1:T % # of frames   
    x = X(find(X(:,t)~=0),t);
    y = Y(find(X(:,t)~=0),t);
    I = INT(find(X(:,t)~=0),t);
    ro = TraceRo(find(X(:,t)~=0),t);
    %z = Z(find(X~=0));
    x = ceil(x);
    y = ceil(y);
    A = single(zeros(yNpx,xNpx));
    A(sub2ind([yNpx,xNpx],x,y)) = ro; % pits
    pitStack2(:,:,t) = A;     
    figure(2)
    imagesc(pitStack2(:,:,t));
end

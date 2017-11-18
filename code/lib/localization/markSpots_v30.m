% puts marks on the spots where they appear in the image stack
    close all;
    isLoadSpots = 1;
    stackOutName = 'detectedSpots_stack.tif';
    load('inputInfo.mat', 'fname')
    d = rdir('traceData-coeff*.mat');
    if isempty(d) || isLoadSpots 
        d = rdir('traceData0-coeff*.mat');
        if isempty(d) || isLoadSpots 
            d = rdir('xyzDataGaus-coeff*.mat');
            if isempty(d) || isLoadSpots 
                d = rdir('xyzData-coeff*.mat');
                if isempty(d) 
                    d = rdir('spotWin-coeff*.mat');
                    if isempty(d)
                        disp('no data file')
                    else
                        load(d.name); lc = 0; % load case
                        disp('loading(0): detected spot center data')
                    end
                else
                    load(d.name); lc = 1; % load case
                    disp('loading(1): center of intensity localized data')
                end
            else
                load(d.name); lc = 2;
                disp('loading(2): Gaus-fit localized data')
                Xg = X; Yg = Y;
                d = rdir('xyzData-coeff*.mat');
                load(d.name); 
            end
        else
            load(d.name); lc = 3;
            disp('loading(3): trace data')
            d = rdir('xyzDataGausFilt-coeff*.mat');
            load(d.name); 
            disp('loading(2): Gaus-fit localized data')
            Xg = X; Yg = Y;
            d = rdir('xyzData-coeff*.mat');
            load(d.name); 
        end
    else
        load(d.name); lc = 4;
        disp('loading(4): combined trace data')
    end
    
    
    
    if lc >=2 % gaus fit
        [trNo frNo] = find(isnan(Xg));
        [trFit frFit] = find(Xg>0);
    elseif lc == 1
        [tr fr] = find(X>0);
    elseif lc == 0
        [tr fr] = find(xBW>0);
        X = xBW; clear xBW;
        Y = yBW; clear yBW;
    end
    if exist(stackOutName), delete(stackOutName); end;
    %if exist('TraceX'), X = TraceX; Y = TraceY; end
    imageInfo=imfinfo(fname);
    Frames=length(imageInfo); imageInfo = imageInfo(1);
    Frames = size(X,2);
    %Frames =10;
    tit = 'title';
    img = imread(fname,1);
    imSz = size(img);
    xx1 = 31; xx2 = xx1+321-1; % crop (sample_1-4 dorsal)
    xx1 = 1; xx2 = imSz(2); % crop (sample_1-4 dorsal)
    yy1 = 1; yy2 = imSz(1); % crop
    mag = 1;
    [magImg, pos, m_, n_] = calcMaxMag(img(yy1:yy2,xx1:xx2),mag);
    dig = floor(log10(Frames))+1;
    
    %X = X*magImg; Y = Y*magImg;
    for ixFrm = 1:Frames
        fig=figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos/2 m_ n_]);
        axe=axes('Parent',fig,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off');
        img = imread(fname,ixFrm);
        if exist('cropCoor.mat')
            load cropCoor; % yy1:yy2,xx1:xx2
        end
        %colormap(gray(65536));
        hImg = imagesc(img(yy1:yy2,xx1:xx2)); %axis image; 
        hold on
        if lc >= 2
            ixTr = find(frNo == ixFrm);
            ixNo = trNo(ixTr); % indice for nonfit spot
            ixTr = find(frFit == ixFrm);
            ixFit = trFit(ixTr); % indice for fit spot (lc=2) or localized spots (lc=1)
            hscat = scatter(X(ixFit,ixFrm),Y(ixFit,ixFrm),'g');
            hscat = scatter(X(ixNo,ixFrm),Y(ixNo,ixFrm),'r');
        else
            ixTr = find(fr == ixFrm);
            ix = tr(ixTr); % indice for nonfit spot
            hscat = scatter(X(ix,ixFrm),Y(ix,ixFrm),10,'*b');
        end
        hold off;
        % print images
        imgFig = getframe(gcf);
        imgOut = imgFig.cdata;
        figPos = get(gcf,'Position');
        
        %imgOut = imgOut((figPos(4)-n)/2+1:(figPos(4)+n)/2,(figPos(3)-m)/2+1:(figPos(3)+m)/2,:);
        
        if ixFrm == 2
            imwrite(imgOut,stackOutName,'Compression', 'none') 
        else
            imwrite(imgOut,stackOutName,'WriteMode','append','Compression', 'none') 
        end
        close(fig)
    end
    return;
    
    %
        set(gcf,'PaperPositionMode','auto')
        resLow = 120; % [dpi] resolution
        rLow = sprintf('-r%i',resLow);
        switch dig
            case 1
                print(sprintf('%s_%01i.tif',frameOutName,ixFrm),'-dtiff',rLow); 
            case 2
                print(sprintf('%s_%02i.tif',frameOutName,ixFrm),'-dtiff',rLow); 
            case 3
                print(sprintf('%s_%03i.tif',frameOutName,ixFrm),'-dtiff',rLow); 
            case 4
                print(sprintf('%s_%04i.tif',frameOutName,ixFrm),'-dtiff',rLow); 
            case 5
                print(sprintf('%s_%05i.tif',frameOutName,ixFrm),'-dtiff',rLow); 
        end
        close(fig)
    tiff2stack(frameOutName,stackOutName)

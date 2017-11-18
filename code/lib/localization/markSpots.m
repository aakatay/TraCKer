% puts marks on the spots where they appear in the image stack
    clear
    close all;
    isLoadSpots = 0;
    isDispOnEveryFrame = 1;
    stackOutName = 'detectedSpots_stack.tif';
    if exist(stackOutName), delete(stackOutName); end;
    
    load fname
    
    %% load images
    imgsFolder = brdir('imgs'); % recursive folder search
    fn = dir('img*.tif');
    imgFile = [];
    if isempty(fn)
        fn = dir([imgsFolder '\img*.tif']);
        if isempty(fn)
            warning('no image file exists')
        else
            imgFile = [imgsFolder '\' fn(1).name];
        end
    else
        imgFile = fn(1).name;
    end
    
    %% load data
    
    d = rdir('traceData-coeff*.mat');
    if isempty(d) || isLoadSpots 
        d = rdir('traceData0-coeff*.mat');
        if 1 || (isempty(d) || isLoadSpots)
            d = rdir('xyzDataGaus-coeff*.mat');
            if isempty(d) || isLoadSpots 
                d = rdir('xyzData-coeff*.mat');
                if isempty(d) || isLoadSpots 
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
    
    %% image info
    tit = 'title';
    img = imread(fname,1);
    sizeImg = size(img);
    Boy1 = sizeImg(1);
    En1 = sizeImg(2);
    binFrame = cfg.img.binFrame;
    Frames = numel(ixSptFrm)-1;
    Frames = floor(Frames/binFrame);
    %dig = floor(log10(Frames))+1;
    %X = X*magImg; Y = Y*magImg;
    hWB =  waitbar(0,'marking spots...');
    frstFrm2 = cfg.img.frstFrm;
    frameVec = 1:Frames-frstFrm2+1;
    %frameVec = 405:412;
    
    %% select XY crop
    isCropXY = 0;
    if isCropXY
        preImg = imread(imgFile);
        if exist('imROI.mat')
            load imROI;
            %xx1 = 30; xx2 = 65; % crop
            %yy1 = 20; yy2 = 65; % crop
        else % select the ROI
            imagesc(preImg);
            hPoly = imrect;
            pos = getPosition(hPoly);
            pos = round(pos);
            posRectCrop = pos;
            xx1 = pos(1);
            yy1 = pos(2);
            xx2 = pos(1)+pos(3)-1; 
            yy2 = pos(2)+pos(4)-1;
            save('imROI','xx1','xx2','yy1','yy2');
            %delete(hPoly)
        end
        En1 = xx2-xx1+1;
        Boy1 = yy2-yy1+1;
    else
        xx1 = 1; xx2 = En1; % crop
        yy1 = 1; yy2 = Boy1; % crop
    end
    mag = 50;
    [magImg0, pos0, m0_, n0_] = calcMaxMag(zeros(Boy1,En1),mag); % recruitment data on image stack
    
    [magImg, ~, ~, ~] = calcMaxMag(zeros(Boy1,En1*2),mag); % recruitment data on pre-image
    [magImg, pos, m_, n_] = calcMaxMag(zeros(Boy1,En1),magImg);
    pos(1) = 1;
    
            %imagesc(preImg(yy1:yy2,xx1:xx2));

    X = X - xx1 + 1;
    Y = Y - yy1 + 1;
    %% read and display image and detections
    for ixFrm = frameVec
        ixSpt = ixSptFrm(ixFrm):ixSptFrm(ixFrm+1)-1; % num spots in the frame;
        for iAv = 1 : binFrame
            frmRead = (ixFrm+frstFrm2-2)*binFrame+iAv;
            temp(:,:,iAv) = imread(fname,frmRead);
        end
        % select the detection in the ROI
        ixSpt2 = ixSpt( find((X(ixSpt)>0) .* (X(ixSpt)<En1) .* (Y(ixSpt)>0) .* (Y(ixSpt)<Boy1) )); %#ok<FNDSB>
        if isDispOnEveryFrame
            % recruitment data on image stack
            fig0=figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos0/2 m0_ n0_]);
            axe=axes('Parent',fig0,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off');

            img = mean(temp,3);
            %colormap(gray(65536));
            %subplot(1,2,1)
            hImg = imagesc(img(yy1:yy2,xx1:xx2)); %axis image; 
            hold on
            if lc >= 2 && ~isempty(ixSpt)
                clr = [ 1  0.5 1];
                clr = [ 1  0 0];
                %hscat = scatter(X(ixSpt),Y(ixSpt),magImg0*3,repmat(clr,numel(ixSpt),1),'.');
                hscat = scatter(X(ixSpt),Y(ixSpt),magImg0*9,repmat(clr,numel(ixSpt),1),'.');
            end
            hold off;
            % print images
            imgFig = getframe(gcf);
            imgOut = imgFig.cdata;
            figPos = get(gcf,'Position');
            if ixFrm == 1
                imwrite(imgOut,stackOutName,'Compression', 'none') 
            else
                imwrite(imgOut,stackOutName,'WriteMode','append','Compression', 'none') 
            end
            close(fig0)
        end
        % display on the pre image
        if 0
            if ixFrm == frameVec(1)
                fig=figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos/2 m_ n_]);
                axe=axes('Parent',fig,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off');
            else
                figure(fig);
            end
            hImg = imagesc(preImg(yy1:yy2,xx1:xx2)); %axis image; 
            hold on
            if lc >= 2 && ~isempty(ixSpt)
                hscat = scatter(X(ixSpt),Y(ixSpt),magImg0,repmat([ 1  0.5 1],numel(ixSpt),1),'*');
            end
            hold off;
        end

        %imgOut = imgOut((figPos(4)-n)/2+1:(figPos(4)+n)/2,(figPos(3)-m)/2+1:(figPos(3)+m)/2,:);

    end
    close(hWB);
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

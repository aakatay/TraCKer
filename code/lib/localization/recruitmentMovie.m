    clear all; close all;

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
    d = rdir('traceJmplessData-coeff*.mat');
    load(d.name);
    d = rdir('traceData-coeff*.mat');
    load(d.name,'cfg');
    
    
    %% image info
    tit = 'title';
    img = imread(fname,1);
    sizeImg = size(img);
    Boy1 = sizeImg(1);
    En1 = sizeImg(2);
    binFrame = cfg.img.binFrame;
    Frames = cfg.img.lastFrm;
    Frames = floor(Frames/binFrame);
    %dig = floor(log10(Frames))+1;
    %X = X*magImg; Y = Y*magImg;
    hWB =  waitbar(0,'marking spots...');
    frstFrm2 = cfg.img.frstFrm;
    frameVec = 1:Frames-frstFrm2+1;
    
    
    %% select XY crop
    preImg = imread(imgFile);
    isCropXY = 1;
    if isCropXY
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
    
    X = X - xx1 + 1;
    Y = Y - yy1 + 1;
    
    %% read and display image and detections
    for ixFrm = frameVec
        for iAv = 1 : binFrame
            frmRead = (ixFrm+frstFrm2-2)*binFrame+iAv;
            temp(:,:,iAv) = imread(fname,frmRead);
        end
        % select the detection in the ROI
        ixSpt2 = ixSpt( find((X(ixSpt)>0) .* (X(ixSpt)<En1) .* (Y(ixSpt)>0) .* (Y(ixSpt)<Boy1) )); %#ok<FNDSB>
        
        % display on the pre image
        if ixFrm == 1
            fig=figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos/2 m_ n_]);
            axe=axes('Parent',fig,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off');
        else
            figure(fig);
        end
        hImg = imagesc(preImg(yy1:yy2,xx1:xx2)); %axis image; 
        hold on
        if lc >= 2 && ~isempty(ixSpt)
            hscat = scatter(X(ixSpt),Y(ixSpt),magImg0,repmat([ 1  0.5 1],numel(ixSpt),1),'.');
        end
        hold off;
    end
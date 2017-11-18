function [intCoef4Lap,smCoeff] = findDiffRec
% calculates the structure replacement
% intMul: illum intensity of laptime acq. wrt bleach acq.

    acqTime=30e-3; % [sec]


    %% find single molecule intensity
    fnameTracerData = rdir('..\_*\traceData0*');
    fnameCoeff = rdir('..\_*\CoeffFit.mat');
    fn = fnameTracerData(1).name;
    %[intSingleMolperFrame1, intSingleMolperTrace ] = findSingleMolIntensity(fn,acqTime);

    load(fnameCoeff.name); % single molecule detection coefficient
    smCoeff = CoeffFit(1);
    
    %% intensity leveling
    intLaser0 = 10; %[au]
    intLaser1 = 100; %[au]
    acqTime0 = 100e-3; %[ms]
    acqTime1 = 30e-3; %[ms]
    intCoef4Lap = intLaser0*acqTime0 / (intLaser1*acqTime1);
    %intSingleMolperFrame0 = intCoef4Lap*intSingleMolperFrame1;

    % FRAPrecoverTime = 100; % [sec]
    % bleachHalfTime = 50; % [sec]
    % acqDuration1 = 90; % [sec]
    % %recoveryRate = int(exp(-acqDuration1/bleachHalfTime)) / FRAPrecoverTime;
    % recoveryRate = 1; % CALCULATE 
    % warning('using recoveryRate = 1, need to be calc.')

    return;
    
    
    
    
    %% load Data
    fnameMapDIR=rdir('structMap*.mat');
    fnameMap = fnameMapDIR(1).name;
    load(fnameMap); % Lt&Bt
    fnameRecDIR = rdir('..\*binImgRcrtSum_time*');
    fnameRec = fnameRecDIR(1).name;

    %% find background
    fnameIllumDIR = rdir('..\*AVGillum*');
    fnameIllum = fnameIllumDIR(1).name;
    imgIllum0 = double(imread(fnameIllum));
    cnvSz=15;
    cnvWin = ones(cnvSz)/cnvSz.^2;
    imgIllum = conv2(imgIllum0,cnvWin,'same');
    imagesc(imgIllum)
    axis image;
    cF = round(cnvSz/2)-1; % convolve frame size
    xSubs=[repmat(imgIllum(:,cF+1),1,cF) repmat(imgIllum(:,end-cF),1,cF)];
    imgIllum(:,[1:cF end-cF+1:end]) = xSubs;
    ySubs=[repmat(imgIllum(cF+1,:),cF,1); repmat(imgIllum(end-cF,:),cF,1)];
    imgIllum([1:cF end-cF+1:end],:) = ySubs;
    imagesc(imgIllum);

    imgIllum=imgIllum/max(imgIllum(:));
    % even illum image
    fnameADIR = rdir('..\*AVG_lap*'); % pre bleach image
    fnameA = fnameADIR(1).name;
    A = double(imread(fnameA));
    A = A./imgIllum; % normalize for uneven illumination
    imagesc(A); axis image

    fnameAN =  ['..\AVG_illumNormLap' fnameA(11:end)]; % normalized A
    fnameAS =  ['..\AVG_strLap' fnameA(11:end)]; % structures in A
    fnameAB =  ['..\AVG_bckgrndLap' fnameA(11:end)]; % structures in A
    imwrite(uint16(A),fnameAN);
    thr = 2100;
    Abw = im2bw(A/thr,1);
    imwrite(uint16(Abw),fnameAS);
    Abckgrnd=A;
    AbwCnv = conv2(double(Abw),ones(3),'same');
    Abckgrnd=Abckgrnd-min(A(:));
    Abckgrnd(AbwCnv>0)=0;
    imwrite(uint16(Abckgrnd),fnameAB);

        %% recruitment image
        imgREC = imread(fnameRec);
        %[gausKernelSz,gausKernelSg] = findSingleMolPSF();
        warning('making up PSF values')
        gausKernelSz = 20;
        gausKernelSg = 4;
        convWin = fspecial('gaussian', gausKernelSz, gausKernelSg);
        imgRECconv0 = conv2(imgREC,convWin,'same');
        % pit image

        imgPits1 = conv2(imgREC,ones(2),'same'); % 2*25nm 50nm
        imgPits2 = conv2(imgREC,ones(4),'same'); % 4*25nm 100nm
        mxRecrt_25nm = max(imgREC(:)); % [#rec/(50nm)^2] 
        mxRecrt_search50nm = max(imgPits1(:)); % [#rec/(50nm)^2] 
        mxRecrt_search100nm = max(imgPits2(:)); % [#rec/(100nm)^2] 

        %imgRECLR = binImage(imgRECconv0,4);
        threshRec = 1; % [rec/16HRpx] THRESHOLD
        convFact = convWin(:);
        convFact = convFact(round(numel(convFact)/2)); % center px value of the convWin
        intThreshPer1HRpx0 = convFact*threshRec;

    %     %% time lapse image
            Ahr = repelem(A,4,4)/16;
            Ahr = Ahr/intSingleMolperFrame0;
    %     stFN = stFNdir(i).name;
    %     if strcmp(stFN(1:6),'STedge'), continue; end
    %     imgSThr = imread(stFN); % high resolution
    %     imgST = binImage(imgSThr,4);
    %     %background
    %     BACKmean=[min(mean(imgSThr,1)),min(mean(imgSThr,2))];
    %     BACK =min(BACKmean);
    %     imgSThr = imgSThr-BACK;
    %     imgSThr(imgSThr<0)=0;

    %%
    for i = 1:numel(Bt)
        close all
        fnameRintDiffHist = dir([ 'rintDiffHist*' sprintf('%02i',i) '.tif' ]);
        if ~isempty(fnameRintDiffHist), continue; end
        lt = Lt(:,:,i);
        [Y,X] = find(lt);
        [y1y2] = minmax(Y'); y1=y1y2(1); y2=y1y2(2);
        [x1x2] = minmax(X'); x1=x1x2(1); x2=x1x2(2);
        sf = 8;
        ry1 = y1 - sf/2;
        rx1 = x1 - sf/2;
        ry2 = y2 + sf/2;
        rx2 = x2 + sf/2;
        y1 = y1 - sf;
        x1 = x1 - sf;
        y2 = y2 + sf;
        x2 = x2 + sf;
        [szY, szX] = size(Ahr);

        if y1<1, y1=1;ry1=1; end; if x1<1, x1=1;rx1=1; end
        if x2>szX, x2=szX;rx2=szX; end; if y2>szY, y2=szY;ry2=szY; end
        sx = rx2-rx1+1;
        sy = ry2-ry1+1;
        a = Ahr(y1:y2,x1:x2);
        imagesc(a)
        r = double(imgREC).*lt;
        [Y,X] = find(r);
        hold on;
        Y = Y - y1+1;
        X = X - x1+1;
        scatter(X,Y,'.','r')
        hr = imrect(gca,[rx1-x1+1 ry1-y1+1 sx sy]);
        hold off

        f = figure(2);set(gcf,'pos',[1230 400 250 150]);
        h = uicontrol('Position', [20 20 100 40], 'String', 'Continue', ...
                      'Callback', 'setRect=1;close(2)');

        figure(1); title('set size of ROI and click continue')
        setRect=0;
        while 1
            pos = getPosition(hr);
            while isempty(find(getPosition(hr)-pos)~=0) && setRect==0
                pause(0.5)
            end
            pause(0.5)

            pos = round(getPosition(hr));
            delete(hr);
            hr=imrect(gca,pos);
            if setRect==1, break; end;
        end
        stWinPos(i,:) = getPosition(hr)+[x1-1 y1-1 0 0];

        imgSThr = imcrop(Ahr,stWinPos(i,:));
        %imgSThr = imgSThr;
        imgRECconv = imcrop(imgRECconv0,stWinPos(i,:));

        %% intensity normalization
        intST=sum(sum(imgSThr));
        intREC=sum(sum(imgRECconv));
        intSumRates = intST/intREC;
        %imgRECconv = uint16(imgRECconv);

        intThreshPer1HRpx = intThreshPer1HRpx0;
        mxRec = double(max([max(imgRECconv(:)) max(imgSThr(:))]));

        % struct name
        posTx = sprintf('%iX%iY%ix%i',pos);
        labelTx = sprintf('_%s_%02i',posTx,i);
        %% THRESHOLD : remove low intensity pixels
            % removal rates
            n1 = numel(find(imgSThr(imgSThr>0)<intThreshPer1HRpx));
            N1 = numel(find(imgSThr>0));
            n2 = numel(find(imgRECconv(imgRECconv>0)<=intThreshPer1HRpx));
            N2 = numel(find(imgRECconv>0));
            r1 = n1/N1*100; % removal rate
            r2 = n2/N2*100;
        ixThreshST = (imgSThr(:)<intThreshPer1HRpx);
        ixThreshREC = (imgRECconv(:)<intThreshPer1HRpx);
        ixThresh = find(ixThreshREC.*ixThreshST);
        imgDiffRec = double(imgRECconv)-double(imgSThr);
        imgSumRec = double(imgRECconv)+double(imgSThr);
        ixNonzero0 = find(imgSumRec(:));
        nPx = numel(ixNonzero0); % total px
        imgSumRec(ixThresh) = 0;
        ixNonzero = find(imgSumRec(:)); % NONZERO PXS
        nTresh = nPx-numel(ixNonzero); % # of removed px
        disp(sprintf('%.02f %% is removed (#%d:%s)',100*nTresh/nPx,i,posTx));

        imgDiffRecLin = imgDiffRec(ixNonzero);
        %% 1/3 histogram 
        imgDiffRecLin = sort(imgDiffRecLin(:));
        imgDiffRecLin = imgDiffRecLin(imgDiffRecLin~=0);
        %figure(5);plot(imgDiffRecLin);
        imgDiffRecLinNorm = imgDiffRecLin;


        figure(6); % histogram
        fnameH = [ 'rintDiffHist' labelTx '.tif'];
        fnameMAT = [ 'rintDiffHist' labelTx '.mat'];
        bw=0.5; % [#recs]binwidth
        [N,xe]=histcounts(imgDiffRecLinNorm,'BinWidth',bw); % xe: edges
        xc = xe(1:end-1)+bw/2; % center pos
        bar(xc,N);
        xcmx = max(abs(xc))+2;
        axis([-xcmx, xcmx, 0, max(N) ])
        grid
        xlabel('intensity difference in # of recruitments');
        ylabel('# of pixels');
        title('intensity difference distribution')
        imgFig = getframe(gcf);
        imgData = imgFig.cdata; 
    %imwrite(uint16(imgData),fnameH)
        matARR = [xc' N'];
        save(fnameMAT,'matARR');

            %% 2/3 overlay images
            figure(7)
            mg=8; % magnification
            %CM=colormap('jet');
            %colormap(mxRec*CM);

            imgCol = [];
            imgCol(:,:,1) = repelem(imgRECconv,mg,mg);
            imgCol(:,:,2) = repelem(imgSThr,mg,mg);
            imgCol(:,:,3) = 0;
            imgCol = imgCol/max(imgCol(:));
            imshow(imgCol);
            szX = size(imgCol,2); szY = size(imgCol,1);
            set(gcf,'units','pixels','Position',[200,200,szX,szY]); 
            set(gca,'units','pixels','Position',[0,0,szX,szY]); 
            imgFig = getframe(gcf);
            imgData = imgFig.cdata; 
            fnameOverlayImg = [ 'imgOverlay' labelTx '.tif' ];
            imwrite(uint16(imgData),fnameOverlayImg);

            %imgCol3 = reshape(imgCol,numel(imgCol)/3,3);
            %sumChan_imgCol3 = sum(imgCol3,1); if sumChan_imgCol3(1)~=sumChan_imgCol3(2), error('normalization error'); end

            fnameOverlayStack = [ 'imgOverlayStack' labelTx '.tif'];
            delete(fnameOverlayStack);
            %% 1: structure
            image(double(imgSThr)/mxRec*64)
            title('structure')
            axis image;
            set(gcf,'units','pixels','Position',[200,200,szX+40,szY+60]); 
            set(gca,'units','pixels','Position',[30,30,szX,szY]); 
            imgFig = getframe(gcf);
            imgData = imgFig.cdata; 
            imwrite(uint16(imgData),fnameOverlayStack,'tiff','WriteMode','append','Compression', 'none');
            %% 2: recruits
            image(double(imgRECconv)/mxRec*64); 
            title('recruitment')
            axis image;
            imgFig = getframe(gcf);
            imgData = imgFig.cdata; 
            imwrite(uint16(imgData),fnameOverlayStack,'tiff','WriteMode','append','Compression', 'none');
            %% 3: difference images1:Static
            diffStruct = -imgDiffRec.*im2bw(-imgDiffRec,0);
            image(diffStruct/mxRec*64); 
            title('STATIC'); % structure-recruitment (>0)
            axis image;
            imgFig = getframe(gcf);
            imgData = imgFig.cdata; 
            imwrite(uint16(imgData),fnameOverlayStack,'tiff','WriteMode','append','Compression', 'none');
            %% 4: difference images2:Growing
            diffRec = imgDiffRec.*im2bw(imgDiffRec,0);
            image(diffRec/mxRec*64)
            title('GROWING'); % recruitment-structure (>0)
            axis image;
            imgFig = getframe(gcf);
            imgData = imgFig.cdata; 
            imwrite(uint16(imgData),fnameOverlayStack,'tiff','WriteMode','append','Compression', 'none');
            %% 5: replacements
            recReplacement = double(imgRECconv)-diffRec;
            image(recReplacement/mxRec*64); 
            title('replacement recruits')
            axis image;
            imgFig = getframe(gcf);
            imgData = imgFig.cdata; 
            imwrite(uint16(imgData),fnameOverlayStack,'tiff','WriteMode','append','Compression', 'none');
            %% 6: pits in center vs edge
            colormap('jet')
            mg = 8;
            imgCol2=[];
            imgCol2(:,:,1) = repelem(recReplacement,mg,mg);
            imgCol2(:,:,2) = repelem(imgSThr,mg,mg);
            imgCol2(:,:,3) = repelem(diffRec,mg,mg);
            imgCol2 = imgCol2/max(imgCol2(:));
            imshow(imgCol2);
            title('pits:white&blue')% static: green, replacement:yellow, pit@edge:blue, pit@center:white
            imgFig = getframe(gcf);
            imgData = imgFig.cdata; 
            imwrite(uint16(imgData),fnameOverlayStack,'tiff','WriteMode','append','Compression', 'none');



        return
        continue;
        % 3/3 linear
        [imgREClin,ix] = sort(imgREC(:));
        imgSTlin = imgST(ix);
        x = 1:numel(imgSTlin);
        figure;plotyy(x,imgREClin,x,imgSTlin)
        cc=3;
        %%
        imgST0 = imgST(:);
        imgST0 = imgST0(imgST0>0);
        numel(imgST0)
        numel(find(imgST0<intSTthreshRecPer16HRpx*1))
        %%
    end

    q=input('do you want to join mat files(y/n)');
    if strcmp(q,'y')
        fnameMATdir = rdir('rintDiffHist*.mat');
        nStr = numel(fnameMATdir);
        szMx =0;
        for i = 1:nStr % check filenames
            fnameMAT = fnameMATdir(i).name;
            labelNo(i) = str2num(fnameMAT(end-5:end-4));
            load(fnameMAT);
            if szMx<size(matARR,1), szMx=size(matARR,1); end
        end
        labelNo = sort(labelNo);
        if ~isempty(find((labelNo-[1:nStr])~=0))
            error('missing or multiple structs');
        end
        matAll = nan(szMx,nStr);
        for i = 1:nStr 
            fnameMAT = fnameMATdir(i).name;
            load(fnameMAT);
            os =(szMx - size(matARR,1))/2; % offset
            matAll(os+1:end-os,:) = matARR(:,2);
        end
        matARR = matAll;
        cellLabel = fnameMap(10:strfind(fnameMap,'.mat')-1);
        matALLfn = sprintf('rintDiffHist%s.mat',cellLabel);
        save(matALLfn,matARR);    
    end

    return
    %%
    figure(18); 
    maximize;
    subplot(1,2,1);
    imagesc(imgPits); 
    colorbar;
    axis image;
    subplot(1,2,2);
    imagesc(imgPits)
    colorbar;
    axis image;

    %% test imread

    colormap('gray')
    [A,cm] = imread(fnameO,'tif');


    %%
            f = @(x) exp(-x.^2).*log(x).^2
            f = @(x) exp(-x)
            Q = integral(f,0,Inf)



    %%
    bleachTime1 = acqTime/intSingleMolperFrame0*intSingleMolperTrace;
    bleachTime0 = bleachTime1*intLaser1/intLaser0;
    collectionRate1 = 1; %  by definition 
    collectionRate0 = acqTime0/bleachTime0; 
end

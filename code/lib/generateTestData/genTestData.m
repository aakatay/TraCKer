function genTestData
    ixRt=1;goData; 
    if ~exist('syntData')
        mkdir('syntData');
    end
    cd('syntData');    
    if ~exist('3D')
        mkdir('3D');
    end
    cd('3D');  
    
    %% param
    addToOldTraces = 0;
    updGenCfg = 1;
    I = 100; %intensity of the spots
    sbPxRt = 9; % has to be odd (sub pixel ratio
    if ~exist('genCfg.mat') || updGenCfg
        %% PSF posint spread function
        spSz = 5; % has to be odd (spotSzie)
        if ~updGenCfg, addToOldTraces = 0;end;
        %PSFmode='cube';
        PSFmode='Gaussian'; sgmXYrat = 0.5; sgmZrat = 0.5;
        imSzIJ = [459,407,45];  % matrix indexing
        nFrame = 50; % # of frames
        nTrace = 5; % # of traces to be defined
        imSzXY = imSzIJ([2 1 3]);  % image indexing
        save('genCfg.mat','spSz','PSFmode','nFrame','nTrace','imSzXY','imSzIJ','','')
    else
        load genCfg.mat;
    end
    if ~rem(sbPxRt*spSz,2)
        display('ERROR : sub pixel magnification ratio (sbPxRt) and spot size (spSz) has to be odd numbers');
        return;
    end
    spSzSubPx = spSz*sbPxRt; %sub pixel spot size
    if strcmp(PSFmode,'cube')
        PSF = ones(spSzSubPx,spSzSubPx,spSzSubPx);
    elseif strcmp(PSFmode,'Gaussian')
        szXY = (spSzSubPx-1)/2;
        szZ = (spSzSubPx-1)/2;
        PSF = Gauss3D(szXY,szZ,sgmXYrat,sgmZrat);
    end
    if ~exist('pval.mat') || (addToOldTraces)
        h=figure;
        maxfig(h,1)
        %% X
        subplot(3,1,1);
        ixTrace = 1;
        while (ixTrace <= nTrace)
            image(zeros(imSzIJ(2),nFrame));
            [fn trace] = ginputExtra; % fn: frame number
            pvalX(ixTrace,:) = dispPoly;
            fnMin(ixTrace) = ceil(min(fn)); if fnMin(ixTrace)<1, fnMin(ixTrace)=1;end;
            fnMax(ixTrace) = floor(max(fn)); if fnMax(ixTrace)>nFrame, fnMax(ixTrace)=nFrame;end;
            if isempty(find(pvalX(ixTrace,fnMin(ixTrace):fnMax(ixTrace)) > imSzXY(1)))
                ixTrace = ixTrace + 1; % all data in the range, progress to next trace acq.
            else
                display('data in undefined region');
            end
        end

        %% Y
        ixTrace = 1;
        while (ixTrace <= nTrace)
            % plot data in other axes5
            subplot(3,1,1);
            image(zeros(imSzIJ(1),nFrame)); hold on;
            plot(pvalX(ixTrace,fnMin:fnMax)); hold off;
            % input positions for the current axis
            subplot(3,1,2);
            image(zeros(imSzIJ(2),nFrame));
            [fn trace] =  ginputExtra;
            pvalY(ixTrace,:) = dispPoly;
            fnMin(ixTrace) = ceil(max([min(fn) fnMin(ixTrace)]));
            fnMax(ixTrace) = floor(min([max(fn) fnMax(ixTrace)]));
            if isempty(find(pvalY(ixTrace,fnMin(ixTrace):fnMax(ixTrace)) > imSzXY(2)))
                ixTrace = ixTrace + 1; % all data in the range, progress to next trace acq.
            else
                display('data in undefined region');
            end
        end
        %% Z
        ixTrace = 1;
        while (ixTrace <= nTrace)
            % plot data in other axes
            subplot(3,1,1);
            image(zeros(imSzIJ(1),nFrame));hold on;
            plot(pvalX(ixTrace,:)); hold off;
            subplot(3,1,2);
            image(zeros(imSzIJ(2),nFrame));hold on;
            plot(pvalY(ixTrace,fnMin:fnMax)); hold off;
            % input positions for the current axis
            subplot(3,1,3);
            image(zeros(imSzIJ(3),nFrame));
            [fn trace] = ginputExtra;
            pvalZ(ixTrace,:) = dispPoly;
            fnMin(ixTrace) = ceil(max([min(fn) fnMin(ixTrace)]));
            fnMax(ixTrace) = floor(min([max(fn) fnMax(ixTrace)]));
            if isempty(find(pvalZ(ixTrace,fnMin(ixTrace):fnMax(ixTrace)) > imSzXY(3))) 
                ixTrace = ixTrace + 1; % all data in the range, progress to next trace acq.
            else
                display('data in undefined region');
            end
        end        
        for ixTrc = 1:nTrace % for each trace define the range of the data
            dataMask = zeros(1,nFrame);
            dataMask(1,fnMin(ixTrc):fnMax(ixTrc)) = 1;
            pvalX(ixTrc,:) = pvalX(ixTrc,:).*dataMask;
            pvalY(ixTrc,:) = pvalY(ixTrc,:).*dataMask;
            pvalZ(ixTrc,:) = pvalZ(ixTrc,:).*dataMask;
        end
        if addToOldTraces
            pvalXn =  pvalX;
            pvalYn =  pvalY;
            pvalZn =  pvalZ;
            load pval.mat;
            pvalX = [pvalX; pvalXn];
            pvalY = [pvalY; pvalYn];
            pvalZ = [pvalZ; pvalZn];
        end
        save('pval.mat', 'pvalX', 'pvalY', 'pvalZ');
    else
        load pval.mat;
    end
    X = pvalX;
    Y = pvalY;
    Z = pvalZ;
    pvalX = round(pvalX); % position of the big window in 3D image
    pvalY = round(pvalY);
    pvalZ = round(pvalZ);
    Xc = X - pvalX; % center position (wrt to center of the big window)
    Yc = Y - pvalY;
    Zc = Z - pvalZ;
    Xc = round(Xc*sbPxRt); % rounded center position in magnified window (wrt to center of the big window)
    Yc = round(Yc*sbPxRt);
    Zc = round(Zc*sbPxRt);
    
    mtrxEmpty = zeros(imSzIJ); % empty matrix
    mtrxEmpty = mtrxEmpty(:);
    imSzXYsubPx = imSzXY*sbPxRt;
    hW = waitbar(0,'registering spot data');
    spts4D =zeros([imSzXY nFrame]);
    for ixFrm = 1:nFrame % for each frame creates a 3D XYZ image stack
        nZeroIx = find(X(:,ixFrm)>0); %nonZeroIx
        for ixTrcNz = 1:length(nZeroIx)
            bW = (spSz+2)*sbPxRt; % size of the big window
            bigWin = zeros(bW,bW,bW); 
            xc = Xc(nZeroIx(ixTrcNz),ixFrm); % shift values of the PSF wrt BigWindow
            yc = Yc(nZeroIx(ixTrcNz),ixFrm);
            zc = Zc(nZeroIx(ixTrcNz),ixFrm);
            bigWin(sbPxRt+1+xc:end-sbPxRt+xc,sbPxRt+1+yc:end-sbPxRt+yc,sbPxRt+1+zc:end-sbPxRt+zc) = PSF*I;
            bigWinDwnSmp = bigWin(1:sbPxRt:end,1:sbPxRt:end,1:sbPxRt:end); % downsample
            spSzHalf = (spSz+1)/2;
            pvalx = pvalX(nZeroIx(ixTrcNz),ixFrm); pvaly = pvalY(nZeroIx(ixTrcNz),ixFrm); pvalz = pvalZ(nZeroIx(ixTrcNz),ixFrm);
            [posPSFx, posPSFy, posPSFz, bigWinDwnSmp] = clipBigWindow;
            spts4D(posPSFx, posPSFy, posPSFz,ixFrm) = bigWinDwnSmp;
        end
        waitbar(ixFrm/nFrame);
    end
    close(hW);
    
    % save stack images at each frame
    hW = waitbar(0,'writing stack images');
    for ixFrm = 1:nFrame
        % max projection
        imgMaxProj(:,:,ixFrm) = max(spts4D(:,:,:,ixFrm),[],3);
        DigitDiff=floor(log10(nFrame));        
        stackName = getStackName;
        
        for ixZ = 1:imSzIJ(3)
            if ixZ == 1
                imwrite(uint16(spts4D(:,:,ixZ,ixFrm)),stackName,'tiff');
            else
                imwrite(uint16(spts4D(:,:,ixZ,ixFrm)),stackName,'tiff','WriteMode','append');
            end
        end
        waitbar(ixFrm/nFrame);
    end
    close(hW);
    stackWrite(imgMaxProj,'maxProjSynt3D.tif')
    
    function [posPSFx, posPSFy, posPSFz, bigWinDwnSmpOut] = clipBigWindow
        % clips the PSF if the spot is on the edge of the image, so that i
        % can be written in 3D image
        posPSFx = pvalx-spSzHalf:pvalx+spSzHalf;
        posPSFy = pvaly-spSzHalf:pvaly+spSzHalf;
        posPSFz = pvalz-spSzHalf:pvalz+spSzHalf;
        nzX = find(posPSFx>0); % non zero indices
        nzY = find(posPSFy>0);
        nzZ = find(posPSFz>0);
        posPSFx = posPSFx(nzX); % non zero indices
        posPSFy = posPSFy(nzY);
        posPSFz = posPSFz(nzZ);
        if size(nzX,2) ~= size(bigWinDwnSmp,1) + size(nzY,2) ~= size(bigWinDwnSmp,2)+ size(nzZ,2) ~= size(bigWinDwnSmp,3)  
            bigWinDwnSmpOut = bigWinDwnSmp(nzX,nzY,nzZ); % clipping
        else
            bigWinDwnSmpOut = bigWinDwnSmp;
        end
    end
        
        
    function stackName = getStackName
            if ixFrm == 1
                DigitDiff=floor(log10(nFrame));
            else
                DigitDiff=floor(log10(nFrame))-floor(log10(ixFrm)); 
            end
            if DigitDiff == 0
            stackName=['stack_' int2str(ixFrm) '.tif'];
            end

            if DigitDiff == 1
                stackName=['stack_0' int2str(ixFrm) '.tif'];
            end

            if DigitDiff == 2
                stackName=['stack_00' int2str(ixFrm) '.tif'];
            end

            if DigitDiff == 3
                stackName=['stack_000' int2str(ixFrm) '.tif'];
            end
    end

    function pval = dispPoly
        P = polyfit(round(fn),trace,3)
        figure;
        pval = polyval(P,1:nFrame);
        plot(pval); hold on;
        plot(fn,trace); hold off
        btn = 0;
        while btn == 0
            btn = waitforbuttonpress;
        end
        close
    end
    

    function [G] = Gauss3D(szXY,szZ,sgmXYrat,sgmZrat)
        %% Normalized 3D Gauss function
        % sgmXY,sgmZ : sigma values [px]
        sgmXY = sgmXYrat*szXY;
        sgmZ = sgmZrat*szZ;
        x = -szXY:szXY; % grid vectors
        z = -szZ:szZ;
        [a b c]= ndgrid(x,x,z); % GRID
        arg   = -(a.*a + b.*b)/(2*sgmXY*sgmXY)-c.*c/(2*sgmZ*sgmZ);
        G     = exp(arg); % exponential
        gMax = max(G(:));
        G = G/gMax;
        %sliceomatic(G)
    end
        

    function [x y] = ginputExtra
    % Author: Lasse Nørfeldt (Norfeldt) 
    % Date: 2012-04-09
        H = gca;
        set(H, 'YLimMode', 'manual'); set(H, 'XLimMode', 'manual');
        set(H, 'YLim', get(H,'YLim')); set(H, 'XLim', get(H,'XLim'));

        xg = []; yg = [];
        i =1;
        while (1)
            [xi yi] = ginput(1);
            xg = [xg xi]; yg = [yg yi];
            if i == 1
                hold on;
                plot(H, xg(i),yg(i),'ro');
            else
                plot(xg([i-1:i]),yg([i-1:i]),'r');
            end
            btn = 0; k = 0;
            while btn == 0
                btn = waitforbuttonpress;
                k = get(gcf,'CurrentCharacter');
            end
            if k == char(13)
                break;
            end
            i = i +1;
        end
        hold off;

        x = xg; y = yg;
    end
end
    close all;
    clear
    F = findall(0,'type','figure','tag','TMWWaitbar'); delete(F)
    isCombDorsalVentral = 0; % run on the top folder (incl. ventral3D and dorsal3D)
    
    if ~isCombDorsalVentral
        inputInfo = dir('inputInfo*.mat');
        load(inputInfo.name);
        dx = PixelSize; 
        dy = dx;
        dz = PlaneDist;
        gausWin = 5; 
        gausSg = 1;

        if ~exist('frameTime')
            frameTime = 0; % seconds 
        end
        isCropXY = 0;
        isCombTraces = 1;
        isObliqueData = 1; % use for oblique PSF (Z position found by XZ localization)
        isWrite2XLS = 0;

        sptJmp = 7; %[pixels] % the max. distance between spots of a trace in consecutive frames (spot jump)
        sptJmpCombination = 5; %2.3;
        sptReAppearTime = 2; %[frames] 
        minTraceLength = 3; % [frames] traces shorter than this value are not plotted

        cfg = struct;
        cfg.trace= struct;
        cfg.trace.sptJmp = sptJmp;
        cfg.trace.minTraceLength = minTraceLength;
        cfg.trace.sptReAppearTime = sptReAppearTime;
        cfg.trace.sptJmpCombination = sptJmpCombination;    

        getDeSkewedPixelSizes = 0;
        if getDeSkewedPixelSizes, DeSkewImage; end; %getDeSkewedPixelSizes (dx, dy, dz)

        BigWindowSize=WindowSize+4;

        %READ the ORIGINAL file
        if exist('fname')
            ;
        elseif exist('fname.mat')
            load fname
            if ~exist(fname)
                fname = sprintf('..\\%s',fname);
            end
        else
            fname = 'data.tif';
            if ~exist(fname), fname2 = rdir(sprintf([fname(1:numel(fname)-4) '*.tif'])); fname=fname2.name; end;
            save('fname','fname')
        end
        imgFrst = imread(fname);
        [Boy1,En1]=size(imgFrst);
        if isCropXY 
            %xx1 = 31; xx2 = xx1+321-1; % crop (sample_1-4 dorsal)
            xx1 = 150; xx2 = xx1+268-1; % crop (sample_1-4 ventral)
            xx1 = 11; xx2 = xx1+356-1; % crop (sample_1-4 dorsal)
            xx1 = 230; xx2 = xx1+94-1; % crop (sample_1-4 dorsal)
            yy1 = 1; yy2 = Boy1; % crop
            yy1 = 71; yy2 = yy1+91; % crop
            En1 = xx2-xx1+1;
            Boy1 = yy2-yy1+1;
            save('cropCoor','xx1','xx2','yy1','yy2');
        elseif exist('crop')
            xx1 = 11; xx2 = xx1+356-1; % crop (sample_1-4 dorsal)
            %xx1 = 150; xx2 = xx1+268-1; % crop (sample_1-4 ventral)
            yy1 = 1; yy2 = Boy1; % crop
            En1 = xx2-xx1+1;
            Boy1 = yy2-yy1+1;
            save('cropCoor','xx1','xx2','yy1','yy2');
        else
            xx1 = 1; xx2 = En1; % crop
            yy1 = 1; yy2 = Boy1; % crop
        end
        [Boy1,En1]=size(imgFrst(yy1:yy2,xx1:xx2));
        imageInfo=imfinfo(fname);
        Frames=length(imageInfo);
Frames = 9; % # of frames to be processed

        coeffMat = dir('coeff*.mat');
        isCoeffLoaded = 0; Coeff = 0;
        if numel(coeffMat)==1
            load(coeffMat.name)
            display('loading coefficient')
            isCoeffLoaded = 1;
        elseif numel(coeffMat) > 1
            display('ERROR : more than one coeff.mat files');
            return;
        end
        assignFileNames

        % resample the max proj to correct PSF shape
        if isObliqueData
            rat = 25/11;
            szX = En1;
            szY = Boy1;
            xx = [1:szX/rat]*rat;
            yy = 1:szY;
            [XX,YY] = meshgrid(xx,yy);
            fnameRS = [fname(1:end-4) '_RS.tif']; % resampled
            if ~exist(fnameRS) % generate resampled maxProj
                for i = 1:numel(imageInfo)
                    imgProj = imread(fname,i);
                    imgProj = double(imgProj(yy1:yy2,xx1:xx2));
                    imgProjRS = interp2(imgProj,XX,YY);
                    imgProjRS = uint16(imgProjRS);
                    imwrite(imgProjRS,fnameRS,'WriteMode','append','Compression', 'none');
                end
                % copyfile(fname,fnameOld);
                % copyfile(fnameRS,fname);

            end
            imgProjRS = imread(fnameRS,1);
            %xx1 = 1; 
            %xx2 = size(imgProjRS,2); % crop
            fnameOrigin = fname;
            fname = fnameRS;
        end

        gaus=fspecial('gaussian', gausWin, gausSg);
        lap=[-1,-1,-1;-1,8,-1;-1,-1,-1];

        %% FIND X & Y
        if ~exist(posDataFileNm) 
            if isObliqueData
                imgPSFfn = dir('PSF*.tif');
                %error('No sample PSF image is found for oblique detection!!!')
                imgPSF = imread(imgPSFfn(1).name);
                imgPSFbckgrnd = imgPSF(1:20,1:20); % crop the PSF background
                bckThrshld = max(imgPSFbckgrnd(:))*1.2; % background threshold
                imgPSF(imgPSF<bckThrshld)=bckThrshld; 
                imgPSF = imgPSF-bckThrshld; % remove background
                imgPSF = imgPSF(10:25,20:40); % crop the PSF data
                filtPSF = double(imgPSF/max(imgPSF(:))); % normalize
                % imagesc(filtPSF)

                % use follow. to set descroption:
                % setTiffDescription('PSF1.tif',tifDescriptConv)
                % tifDescriptConv = sprintf('bckground: x=1:20, y=01:20 \nPSF: x=10:25, y=20:40');


                % read rotation angle
                fAngle = fopen('rotAngle.txt','r'); 
                rotAngle = fread(fAngle);
                rotAngle = str2num(char(rotAngle'));
                fclose(fAngle);
            end

            if ~isCoeffLoaded
                tit = 'determine the coefficient for background filtering';
                hFig=figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256));
                set(gcf,'units','normalized','outerposition',[0 0 1 1])
                kk =0.1;
                axe=axes('Parent',hFig,'DataAspectRatio',[1 1 1],'Position',[kk kk 1-2*kk 1-2*kk ],'Visible','off','outerposition',[0 0 1 1]);

                % TileParameters
                mag = 2; % magnification
                temp=get(0,'ScreenSize');
                screenX = temp(3); screenY = temp(4);
                nX = round((screenX + 200)/En1/mag);
                nY = round((screenY - 0)/Boy1/mag);
                nFrTile = nX*nY; % number of frames to be tiled
                nRow = nY;  

                for j=1:nFrTile
                    temp=imread(fname,j);
                    if ~isObliqueData % cropping is not needed in oblique
                        PREdata(:,:,j) = temp(yy1:yy2,xx1:xx2);
                    else
                        PREdata(:,:,j) = temp;
                    end
                    PREdataFilt = imfilter(double(PREdata),gaus,'symmetric');
                    PREdataFilt = imfilter(PREdataFilt,lap,'symmetric');
                end
                PREdataTiled = tileFrames(PREdata ,nRow);
                PREdataFiltTiled = tileFrames(PREdataFilt ,nRow);

                SHOW=PREdataTiled;
                imgMx = max(PREdataFiltTiled(:));
                Coeff = imgMx/5;

                hTextCoeff1 = uicontrol('style','text','BackgroundColor',[1 1 1],'String',sprintf('Coefficient: %i',round(Coeff)),'Position',[20 150 160 15]);
                hTextCoeff2 = uicontrol('style','text','BackgroundColor',[1 1 1],'String',{'select the coefficient'; 'use keys to change threshold value';'w: ++0.1 e:--0.1';'s: ++0.01 d:--0.01';'press q to continue'},'Position',[20 20 160 115]);
                q=0; 
                slideVal = double(double(Coeff)/double(imgMx));
                while q == 0
                    colormap(pink);
                    imagesc(PREdataTiled);axis image;
                    PREdataFiltTiledDiv=PREdataFiltTiled/Coeff;
                    BINAR=im2bw(PREdataFiltTiledDiv,1);
                    DINAR=uint16(BINAR).*PREdataTiled;
                    BWSHOW=imregionalmax(DINAR, 8);
                    [y,x,v]=find(BWSHOW==1);
                    figure(1);
                    %imagesc(BWSHOW); % filtered spots
                    hold on; scatter(x,y,22,'o');
                    hold off;
                    set(hTextCoeff1,'String',sprintf('Coefficient: %i',round(Coeff)));
                    axis image;
                    btn = 0;
                    while btn == 0
                        btn = waitforbuttonpress;
                        k = get(hFig,'CurrentCharacter');
                    end
                    upd= 0;
                    switch lower(k)
                        case 'x'
                            slideVal = slideVal - 0.001; upd = 1;
                        case 'c'
                            slideVal = slideVal + 0.001; upd = 1;
                        case 's'
                            slideVal = slideVal - 0.01; upd = 1;
                        case 'd'
                            slideVal = slideVal + 0.01; upd = 1;
                        case 'w'
                            slideVal = slideVal - 0.1; upd = 1;
                        case 'e'
                            slideVal = slideVal + 0.1; upd = 1;
                        case 'q'
                            q = 1;
                    end
                     if upd == 1
                         Coeff = slideVal*imgMx;
                     end
                end
                save('Coeff.mat','Coeff');
                assignFileNames
            end
            close;

            for j=1:Frames
                img =imread(fname,j);
                if ~isObliqueData % cropping is not needed in oblique
                    J(:,:,j) = img(yy1:yy2,xx1:xx2);
                    img = img(yy1:yy2,xx1:xx2);
                else
                    J(:,:,j) = img;
                    imgOrigin =imread(fnameOrigin,j);
                    Jorigin(:,:,j) = imgOrigin;
                end
                img = double(img);

                JF(:,:,j) = imfilter(img,gaus,'symmetric');
                JF(:,:,j) = imfilter(JF(:,:,j),lap,'symmetric');
            end
            J=uint16(J);
            h = waitbar(0,'3D localization...');
            %% find 3D intensity
            sp = 0;
            nSp = 0;
            Coeff = double(Coeff);
            for k=1:Frames
                if exist('BREAK.MAT')
                    if exist('breakVal.mat')
                        clear
                        load breakVal;
                        delete('BREAK.MAT','breakVal')
                    else
                        save('breakVal');
                        return;
                    end
                end
                IMG = J(:,:,k);
                % find peaks
                dataFilt = JF(:,:,k); 
                dataFiltDiv = dataFilt/Coeff;
                bin = im2bw(dataFiltDiv,1);
                if isempty(find( bin == 1)),continue;end
                din = uint16(bin).*IMG;
                BW = imregionalmax(din, 8);
                [B,L] = bwboundaries(BW,'noholes');
                if isObliqueData % cropping is not needed in oblique
                   IMG = Jorigin(yy1:yy2,xx1:xx2,k);
                end

                % remove spots at the boundaries
                if exist('boundaries')
                    scaleOblique=1;
                    if isObliqueData, scaleOblique=211/481;  end;% should be ~1/(size(imgSpotRot,1)/size(imgSpot,1));
                    boundariesScaled = boundaries.*[scaleOblique 1 scaleOblique 1];
                    xb1 = boundariesScaled(1)+1;
                    xb2 = boundariesScaled(1)+boundariesScaled(3);
                    numSpots = numel(B);
                    i = 1;
                    while i<numSpots
                        c=cell2mat(B(i));
                        Py=uint16(mean(c(:,1)));
                        Px=uint16(mean(c(:,2)));
                        if (Px <  xb1) ||  (Px > xb2) % if at the boundary
                            if i == 1
                                B = B([2:end]);
                            else
                                B = B([1:i-1 i+1:end]); % neglect the spot 
                            end
                            i = i - 1;
                            numSpots = numSpots - 1;
                        end
                        i = i + 1;
                    end
                end

                %% read stack
                if StackNum > 1
                    clear Stack;
                    for l=1:StackNum
                        DigitDiff=floor(log10(Frames))-floor(log10(k));   
                        
                        if k == 1
                        DigitDiff=floor(log10(Frames));
                        end
DigitDiff = 2;
                        if DigitDiff == 0   
                        Stack_=imread(['stack_' int2str(k) '.tif'],l);
                        end

                        if DigitDiff == 1   
                        Stack_=imread(['stack_0' int2str(k) '.tif'],l);
                        end

                        if DigitDiff == 2   
                        Stack_=imread(['stack_00' int2str(k) '.tif'],l);
                        end

                        if DigitDiff == 3   
                        Stack_=imread(['stack_000' int2str(k) '.tif'],l);
                        end

                        % Stack_=Stack_-min(min(Stack_));

                        Stack(:,:,l) = Stack_(yy1:yy2,xx1:xx2);
                    end
                    Stack=double(Stack);
                    if isObliqueData, Stack = permute(Stack,[3 2 1]);  end; % YXZ -> ZXY
                end

                q=0; 
                %DEFINE Size
                %[En,Boy]=size(IMG);
                nSp = nSp + length(B);
                for m=1:length(B) % for each spot
                    c=cell2mat(B(m));
                    %csize=(max(c(:,1))-min(c(:,1)))*(max(c(:,2))-min(c(:,2)));
                    q=q+1;sp = sp + 1;
                    Py=uint16(mean(c(:,1)));
                    Px=uint16(mean(c(:,2)));
                    if isObliqueData
                        Px = Px*rat;
                    end
                    isDebug = 0;
                    if isDebug
                        imagesc(IMG);
                        hold on
                        scatter(Px,Py)
                        hold off
                    end
                    % adjust the big window center position for the spots at
                    % the edges
                    PxBW = Px;
                    PyBW = Py;

                    % window center
                    if (Px-(WindowSize+1)/2)<1
                        Px=(WindowSize+1)/2;
                    end
                    if (Py-(WindowSize+1)/2)<1
                        Py=(WindowSize+1)/2;
                    end
                    if (Px+(WindowSize+1)/2)>En1
                        Px=En1-(WindowSize+1)/2;
                    end
                    if (Py+(WindowSize+1)/2)>Boy1
                        Py=Boy1-(WindowSize+1)/2;
                    end
                    % big window center
                    if (PxBW-(BigWindowSize+1)/2)<1
                        PxBW=(BigWindowSize+1)/2;
                    end
                    if (PyBW-(BigWindowSize+1)/2)<1
                        PyBW=(BigWindowSize+1)/2;
                    end
                    if (PxBW+(BigWindowSize+1)/2)>En1
                        PxBW=En1-(BigWindowSize+1)/2;
                    end
                    if (PyBW+(BigWindowSize+1)/2)>Boy1
                        PyBW=Boy1-(BigWindowSize+1)/2;
                    end
                    %DEFINE Window
                    yy = Py-(WindowSize+1)/2;
                    y1 = yy + 1; y2 = yy + WindowSize;
                    xx = Px-(WindowSize+1)/2;
                    x1 = xx + 1; x2 = xx + WindowSize;
                    Window=IMG(y1:y2,x1:x2);

                    %DEFINE Big Window
                    yy = PyBW-(BigWindowSize+1)/2;
                    y1 = yy + 1; y2 = yy + BigWindowSize;
                    xx = PxBW-(BigWindowSize+1)/2;
                    x1 = xx + 1; x2 = xx + BigWindowSize;
                    BigWindow=IMG(y1:y2,x1:x2);

                    %background
                    BACKmean=[min(mean(BigWindow,1)),min(mean(BigWindow,2))];
                    BACK(q,k)=min(BACKmean);

                    %FIND Total Intensity
                    INT(q,k)=sum(sum(Window))-BACK(q,k)*(WindowSize)^2;

                    %FIND Intensity Center
                    TopX=0;
                    TopY=0;
                    TopColum=0;
                    TopRow=0;
                    WSumX=0;
                    WSumY=0;

                    for j=1:WindowSize
                       TopX(j)=sum(Window(:,j));
                    end
                    TopX=TopX-min(TopX);
                    TopRow=sum(TopX);

                    for j=1:WindowSize
                        WSumX=WSumX+j*TopX(j);
                    end

                    for i=1:WindowSize
                       TopY(i)=sum(Window(i,:));
                    end
                    TopY=TopY-min(TopY);
                    TopColum=sum(TopY);

                    for i=1:WindowSize
                        WSumY=WSumY+i*TopY(i);
                    end

                    Xc(k)=WSumX/TopRow;
                    Yc(k)=WSumY/TopColum;
                    if isnan(Xc(k)), Xc(k)=(WindowSize+1)/2; end;
                    if isnan(Yc(k)), Yc(k)=(WindowSize+1)/2; end;


                    PXc=uint8(Xc(k));
                    PYc=uint8(Yc(k));

                    %center of intensity
                    X_=double(Px)+Xc(k)-double((WindowSize+1)/2);
                    Y_=double(Py)+Yc(k)-double((WindowSize+1)/2);

                    X(q,k)=X_;
                    Y(q,k)=Y_;
                    %Inten(Py,Px)=INT(k);
                    if X(q,k)<0
                        pause;
                    end

                    if StackNum > 1
                        %% FIND Z NOW
                        if isObliqueData
                            isdebug = 0;
                            Xyz = round(X(q,k));
                            xYz = round(Y(q,k));
                            dxCrop = round(size(filtPSF,2)*1/4);
                            dzCrop = round(size(filtPSF,1)*1/3);
                            xMin = Xyz-dxCrop;
                            xMax = Xyz+dxCrop;
                            yMin = xYz-(WindowSize-1)/2;
                            yMax = xYz+(WindowSize-1)/2;
                            %yMin = xYz;
                            %yMax = xYz;
                            if xMin < 1
                                xMin = 1;
                            end
                            if xMax > En1
                                xMax = En1;
                            end
                            if yMin < 1
                                yMin = 1;
                            end
                            if yMax > Boy1
                                yMax = Boy1;
                            end
                            imgSpotXZY = Stack(:,xMin:xMax,yMin:yMax);  % Stack in ZXY
                            %figure(1);imagesc(imgSpotXZY)
                            imgSpotXZ = max(imgSpotXZY,[],3);
                            % imagesc(imgSpotXZ)
                            [~, ixZ] = max(max(imgSpotXZ,[],2)); % z coordinate
                            if isdebug
                                figure(20); maximize
                                subplot(4,5,1)
                                imagesc(max(Stack,[],3)); title('stackXZproj'); axis equal; axis tight
                                subplot(4,5,2)
                                imagesc(max(imgSpotXZY,[],3)); title('XY crop'); axis equal; axis tight
                                hold on; 
                                scatter(dxCrop+1,ixZ);
                                hold off;
                            end

                            % localize in Z
                            yukseklik = size(imgSpotXZ,1);
                            zMin = ixZ-dzCrop;
                            zMax = ixZ+dzCrop;
                            imgSpot = zeros(dzCrop*2,size(imgSpotXZY,2),size(imgSpotXZY,3));
                            if zMin < 1
                                zMin = 1;
                            end
                            if zMax > yukseklik 
                                %dzMax = zMax-yukseklik;
                                zMax = yukseklik;
                            end
                            padSize = dzCrop*2-(zMax-zMin)-1;
                            if padSize<0, padSize=0; end
                            imgSpot_ = imgSpotXZY(zMin:zMax,:,:);  % Stack in ZXY
                            imgSpot([zMin:zMax]-zMin+1,:,:) = imgSpot_;

                            %figure(2);imagesc(imgSpot)
                            %[m k]
                            %size(imgSpot)
                            clear imgSpotRot;
                            for i = 1:size(imgSpot,3) 
                                imgSpotRot(:,:,i) = imrotate(imgSpot(:,:,i),rotAngle,'bicubic'); %
                            end

                            if isdebug
                                subplot(4,5,3)
                                imagesc(max(imgSpot,[],3)); title('ZXY crop'); axis equal; axis tight
                                subplot(4,5,4)
                                imagesc(max(imgSpotRot,[],3)); title('rotated'); axis equal; axis tight
                            end

                            [szRtZ, szRtX,~] = size(imgSpotRot);
                            cRtZ = round(szRtZ/2); 
                            cRtX = round(szRtX/2);
                            dxCr = 3;
                            dzCr = 5;
                            ZWindow = imgSpotRot(cRtZ-dzCr:cRtZ+dzCr,cRtX-dxCr:cRtX+dxCr,:);
                            if isdebug
                                subplot(4,5,5)
                                imagesc(max(ZWindow,[],3)); title('ZX crop'); axis equal; axis tight
                            end
                            ZWindow2 = ZWindow;
                            threshZ = max(ZWindow2(:))/1.2;
                            ZWindow2(ZWindow2<threshZ)=threshZ;
                            CrSzZ = dzCr*2+1; % crop size Z

                            % same localization procedure
                            ZWindowXYZ = permute(ZWindow2,[2 3 1]); % --> XYZ
                            TopZr=squeeze( sum(sum(ZWindowXYZ,2),1)-sum(min(ZWindowXYZ)*WindowSize,2));
                            zr = 1:CrSzZ;
                            TopZr=TopZr-min(TopZr); % remove background
                            Z_=CrSzZ-sum(zr.*TopZr')/sum(TopZr);
                            %Z(q,k)=Z_+cRtZ-dzCr/2;
                            rotRat = sin(rotAngle/180*pi); % should be ~1/(size(imgSpotRot,1)/size(imgSpot,1));
                            Z(q,k) = size(imgSpotXZY,1)-((cRtZ-dzCr+Z_-1)*rotRat-padSize+yukseklik-zMax);

                            if isdebug
                                subplot(4,5,6)
                                imagesc(max(ZWindow2,[],3)); title('ZX crop (even background)'); axis equal; axis tight
                                hold on
                                scatter(dxCr+1,CrSzZ-Z_);
                                hold off

                                subplot(4,5,9)
                                imagesc(max(imgSpotRot,[],3)); title('rotated'); axis equal; axis tight
                                hold on
                                scatter(round(size(imgSpotRot,2)/2),size(imgSpotRot,1)-(cRtZ-dzCr+Z_-1));
                                hold off


                                subplot(4,5,8)
                                imagesc(max(imgSpot,[],3)); title('ZXY crop'); axis equal; axis tight
                                hold on
                                scatter(round(size(imgSpot,2)/2),size(imgSpot,1)-(cRtZ-dzCr+Z_-1)/rotRat);
                                hold off

                                %
                                subplot(4,5,7)
                                imagesc(max(imgSpotXZY,[],3)); title('XY crop'); axis equal; axis tight
                                hold on
                                scatter(round(size(imgSpotXZY,2)/2),Z(q,k));
                                hold off

                                figLabel = sprintf('frame%i, spot%i',k,m);
                                set(gcf,'Name',figLabel)

                                imgFig = getframe(gcf);
                                imgOut = imgFig.cdata;
                                imwrite(imgOut,[figLabel '.tif'],'Compression', 'none') 

                                subplot(4,5,10)
                                imagesc(BigWindow); title('stackXYproj'); axis equal; axis tight
                                imagesc(IMG); title('stackXYproj'); axis equal; axis tight
                                hold on
                                scatter(X(q,k),Y(q,k));
                                hold off
                                dd = 20;
                                set(gca,'Xlim',[X(q,k)-dd X(q,k)+dd],'Ylim',[Y(q,k)-dd Y(q,k)+dd])

                                subplot(4,5,[11:15])
                                imagesc(max(Stack,[],3)); title('stackXZproj'); axis equal; axis tight
                                hold on
                                scatter(X(q,k),Z(q,k));
                                hold off

                                subplot(4,5,[11:15]+5)
                                imagesc(Stack(:,:,xYz)); title('stackXZproj'); axis equal; axis tight
                                hold on
                                scatter(X(q,k),Z(q,k));
                                hold off
                            end
                        else
                            ZWindow = Stack(Py-(WindowSize+1)/2+1:Py+(WindowSize+1)/2-1,Px-(WindowSize+1)/2+1:Px+(WindowSize+1)/2-1,:);
                            TopZr=squeeze( sum(sum(ZWindow,2),1)-sum(min(ZWindow)*WindowSize,2));
                            % older ver. : TopZr(l)=sum(sum(ZWindow))-sum(min(ZWindow))*WindowSize; %  data range in the XY plane
                            zr = 1:StackNum;

                            z=zr;
                            TopZ=TopZr';
                            TopZ=TopZ-min(TopZ); % remove background
                            % fz=fit(z,TopZ,'gauss1');
                            % coeffvalues(fz);
                            % Z(q,k)=ans(2);
                            Z(q,k)=StackNum-sum(z.*TopZ)/sum(TopZ);

                            zPosZero = 0;
                            if TopZ==0
                                %fprintf('z position is zero check the code');
                                zPosZero = zPosZero + 1;
                                Z(q,k)=StackNum/2;
                                %return;
                            end
                        end
                    else
                        Z = zeros(size(X));
                    end
                end
                waitbar(k/Frames)
                if ~rem((k-1),100) 
                    time100frame = 0;
                    if k > 1
                        time100frame = toc; 
                        disp(sprintf('time for 100 frames : %.02f\n',time100frame));
                    end;
                    tic;

                end
                isDebug1 = 0;
                if isDebug1 % localization
                    imagesc(IMG);
                    hold on
                    scatter(X(:,k),Y(:,k));
                    hold off
                end

                isDebug2 = 0;
                if isDebug2 % localization
                    imagesc(max(Stack,[],3));
                    hold on
                    scatter(X(:,k),Z(:,k));
                    hold off
                end            
            end
            %stackWrite(NBINs*3000,'NBINs.tif')
            clear JF J ;
            close(h);
            save(xyzDataFileNm,'X','Y','Z','INT'); 
            clear X Y Z INT spotWin NBINs Zstack maxProj;
            clear NBIN NDIN BW B L IMG NBINsum temp imgFrst Frames Frames2 ZWindow Stack
            save(posDataFileNm); % save data  


            %% move files
            folderNm= sprintf('%s-coeff%i',label,round(Coeff));
            if exist(folderNm)
                fn = 1;
                while exist(sprintf('%s-coeff%i_%i',label,round(Coeff),fn))
                    fn = fn + 1;
                end
                folderNm = sprintf('%s-coeff%i_%i',label,round(Coeff),fn);
                mkdir(folderNm);
            else
                mkdir(folderNm);
            end

            % move files
            isCopy = 1;
            if isCopy
                copyfile('Coeff.mat',strcat(folderNm,'\',nmCoeff));
                copyfile(fname,strcat(folderNm,'\',fname));      
                if exist('fnameOrigin'), copyfile(fnameOrigin,strcat(folderNm,'\',fnameOrigin)); end
                if exist('cropCoor.mat'), movefile('cropCoor.mat',strcat(folderNm,'\cropCoor.mat')); end;
                movefile(posDataFileNm,strcat(folderNm,'\',posDataFileNm));
                movefile(xyzDataFileNm,strcat(folderNm,'\',xyzDataFileNm));
                if exist('fname.mat')
                    copyfile('fname.mat',strcat(folderNm,'\fname.mat'));  
                else
                    copyfile('inputInfo.mat',strcat(folderNm,'\inputInfo.mat'));         
                end
                delete(fnameRS);
                cd(folderNm)
            end
            %return;
        end

        %% TIME To FIND OUT THE TRACES
        if ~exist(traceDataFileNm0) 
            load(xyzDataFileNm);
            load(posDataFileNm, 'cfg')

            cfg.trace= struct;
            cfg.trace.sptJmp = sptJmp;
            cfg.trace.minTraceLength = minTraceLength;
            cfg.trace.sptReAppearTime = sptReAppearTime;
            cfg.trace.sptJmpCombination = sptJmpCombination;   


            isCropFrames = 0;
            if isCropFrames
                f1=1;f2=300;
                X = X(:,f1:f2);
                Y = Y(:,f1:f2);
                Z = Z(:,f1:f2);
            end
            Xilk=X;
            Yilk=Y;
            Zilk=Z;
            [Boy,Frames]=size(X);

            h = waitbar(0,'Finding the traces...');
            p=0;
            f=0;
            tic
            isDebug = 0;
            for k=1:Frames-1 % number of frames
                if isDebug, disp(sprintf('=frame#:%i\n',m)); end
                for m=1:Boy % number of spots
                    if X(m,k)>0
                    if X(m,k)<En1
                        tracex=zeros(1,Frames);
                        tracey=zeros(1,Frames);
                        tracez=zeros(1,Frames);
                        traceint=zeros(1,Frames);

                        dif=Inf(Boy,Frames-k+1);
                        difbin=zeros(Boy,Frames-k+1);            
                        AslX=X(m,k);       % last X value in the tracking
                        AslY=Y(m,k);       % last Y value in the tracking
                        AslZ=Z(m,k);
                        if isDebug, disp(sprintf('==spot#:%i\n',m)); end
                        quit = 0;
                        ll=1;
                        for l=1:Frames-k+1 % later frames                        
                            if isDebug, disp(sprintf('===check frame#:%i\n',l)); end
                            if quit, 
                                break; 
                            end
                            % distance betw. all spots
                            dif(:,l)=sqrt((AslX-X(:,k+l-1)).^2*dx^2 + (AslY-Y(:,k+l-1)).^2*dy^2+(AslZ-Z(:,k+l-1)).^2*dz^2*(PlaneDist/PixelSize)^2);            

                            [v n2] = min(dif(:,l));
                            BOY = Boy;
                            BOY = 1;
                            for n=1:BOY % all spots
                                n=n2;
                                if l-ll>sptReAppearTime % if the next frame where a spots re-appears in sptJmp distance is 4 frames apart ignores it.
                                    %if sum(sum(difbin(:,l-2:l-1))) == 0, 
                                        quit =1;
                                        break, 
                                    %end
                                end

                                if isDebug, disp(sprintf('===check spot#:%i\n',n)); end
                                if sqrt((AslX-X(n2,k+l-1))^2+(AslY-Y(n2,k+l-1))^2+[(AslZ-Z(n2,k+l-1))*PlaneDist/PixelSize]^2) < sptJmp;
                                    if dif(n,l)==min(dif(:,l))
                                        if n ~= n2
                                            disp(sprintf('n:%i, n2"%i\n',n,n2));
                                            %disp('strange')
                                            break
                                        end         
                                        ll=l;
                                        difbin(n,l)=1;
                                        AslX=Xilk(n,k+l-1);
                                        AslY=Yilk(n,k+l-1);
                                        AslZ=Zilk(n,k+l-1);
                                        break
                                    end
                                end
                            end % spots
                        end % frames

                        for n=1:Boy
                        for l=1:Frames-k+1
                            if difbin(n,l)==1;
                                tracex(k+l-1)=Xilk(n,k+l-1);
                                tracey(k+l-1)=Yilk(n,k+l-1);
                                tracez(k+l-1)=Zilk(n,k+l-1);
                                traceint(k+l-1)=INT(n,k+l-1);
                                X(n,k+l-1)=Inf;
                                Y(n,k+l-1)=Inf;
                                Z(n,k+l-1)=Inf;
                            end
                        end % frames
                        end % spots

                        num=numel(find(tracex>0));
                        if num>=minTraceLength % if # of data points larger than minTraceLength, than saves as a trace
                             pos=find(tracex>0);
                             ilk=zeros(1,num+1);
                             son=zeros(1,num+1);
                             ilk(1:1:num)=pos(1:1:num);
                             son(2:1:num+1)=pos(1:1:num);
                             fark=ilk-son;
                             %if numel(find(fark==1))>2 % # of consecutive data points
                            p=p+1;
                            TraceX(p,:)= sparse(tracex);
                            TraceY(p,:)=sparse(tracey);
                            TraceZ(p,:)=sparse(tracez);
                            TraceINT(p,:)=sparse(traceint);   


                             %end
                        end
                    end % x<EN
                    end % x>0
                end
               waitbar(k / Frames)
               if exist('_stopRunning-ON')
                   break
               end
            end
            TrackTime = toc;
            save('TrackTime','TrackTime');
            close(h)
            %return;
            save(traceDataFileNm0,'TraceX','TraceY','TraceZ','TraceINT','cfg'); 
            delete(traceDataFileNm)
        end


        if ~exist(traceDataFileNm) && isCombTraces   % speed and trace combination
            %% COMBINE TRACES
            load(traceDataFileNm0)
            [m n]=size(TraceX);
            h = waitbar(0,'Combining the traces...');
            disp('Combining the traces...');

            [Boy2,Frames]=size(TraceX);
            tic;
            numCombTraces=0;
            isDebugTrComb = 0;
            for i=1:Boy2 % each trace
            %     TraceDif1=0;TraceDif2=0;TraceDif3=0;
                if i <= Boy2-1  % combine traces
                    LastElement=max(find(TraceX(i,:)>0));
                    LastBefore=LastElement-1;
                    for j=i+1:Boy2 
                        FirstElement=min(find(TraceX(j,:)>0));
                        FirstAfter=FirstElement+1;
                        if FirstElement-LastElement>-1
                            if FirstElement-LastElement<3
                                TraceDif1=sqrt([TraceX(i,LastElement)-TraceX(j,FirstElement)]^2+[TraceY(i,LastElement)-TraceY(j,FirstElement)]^2+[(TraceZ(i,LastElement)-TraceZ(j,FirstElement))*PlaneDist/PixelSize]^2);
                                TraceDif2=sqrt([2*TraceX(i,LastElement)-TraceX(i,LastBefore)-TraceX(j,FirstElement)]^2+[2*TraceY(i,LastElement)-TraceY(i,LastBefore)-TraceY(j,FirstElement)]^2+[(2*TraceZ(i,LastElement)-TraceZ(i,LastBefore)-TraceZ(j,FirstElement))*PlaneDist/PixelSize]^2);
                                TraceDif3=sqrt([TraceX(i,LastElement)-2*TraceX(j,FirstElement)+TraceX(j,FirstAfter)]^2+[TraceY(i,LastElement)-2*TraceY(j,FirstElement)+TraceY(j,FirstAfter)]^2+[(TraceZ(i,LastElement)-2*TraceZ(j,FirstElement)+TraceZ(j,FirstAfter))*PlaneDist/PixelSize]^2);
                                TraceDif=[TraceDif1,TraceDif2,TraceDif3];
                                %TraceDif=sqrt((TraceX(i,LastElement)-TraceX(j,FirstElement))^2+(TraceY(i,LastElement)-TraceY(j,FirstElement))^2);
                                if min(TraceDif)<sptJmpCombination
                                    if isDebugTrComb, scatter(TraceX(i,TraceX(i,:)>0),TraceY(i,TraceY(i,:)>0),'*r');hold on; end;
                                    for k=1:Frames
                                        if [TraceINT(i,k)+TraceINT(j,k)] > 0
                                        TraceX(i,k)= [TraceX(i,k)*TraceINT(i,k)+TraceX(j,k)*TraceINT(j,k)]/[TraceINT(i,k)+TraceINT(j,k)];
                                        TraceY(i,k)= [TraceY(i,k)*TraceINT(i,k)+TraceY(j,k)*TraceINT(j,k)]/[TraceINT(i,k)+TraceINT(j,k)];
                                        TraceZ(i,k)= [TraceZ(i,k)*TraceINT(i,k)+TraceZ(j,k)*TraceINT(j,k)]/[TraceINT(i,k)+TraceINT(j,k)];
                                        TraceINT(i,k)= [TraceINT(i,k)*TraceINT(i,k)+TraceINT(j,k)*TraceINT(j,k)]/[TraceINT(i,k)+TraceINT(j,k)];
                            %             TnaceX(i,k)=TraceX(i,k); TnaceY(i,k)=TraceY(i,k); TnaceX(j,k)=TraceX(j,k); TraceY(j,k)=TraceY(j,k);
                                        end
                                    end
                                    if isDebugTrComb, 
                                        scatter(TraceX(j,TraceX(j,:)>0),TraceY(j,TraceY(j,:)>0),'+b')
                                        scatter(TraceX(i,TraceX(i,:)>0),TraceY(i,TraceY(i,:)>0),'.g')  
                                        hold off
                                    end
                                    TraceX(j,:)=NaN; TraceY(j,:)=NaN; TraceZ(j,:)=NaN; TraceINT(j,:)=NaN;
                                    LastElement=max(find(TraceX(i,:)>0));
                                    LastBefore=LastElement-1;
                                    numCombTraces=numCombTraces+1;
                                end
                            end
                        end
                    end
                end % i <= Boy2-1
                    %Tcomp0 = toc;
                waitbar(i / Boy2)
            end



            close(h)
            % 'traceData-coeff%d.mat'
            Tcomp = toc;
            save(traceDataFileNm,'TraceX','TraceY','TraceZ','TraceINT','numCombTraces','cfg')
            delete(traceDataFileNm2)
        end

        if ~exist(traceDataFileNm2)    % speed and trace combination
        %%  SPEED : fill gaps (missing frames in jumpy traces)    
            if isCombTraces
                load(traceDataFileNm);
            else
                load(traceDataFileNm0);
                numCombTraces = 0;
            end
            [Boy2,Frames]=size(TraceX);
            for i=1:Boy2 % each trace
                tracex = TraceX(i,:);
                tracey = TraceY(i,:);
                tracez = TraceZ(i,:);
                if isnan(tracex(1)), continue; end;

                ind = find(tracex~=0);
                frst = min(ind); last = max(ind);
                jmp = abs((tracex(frst:last)>0)-1);
                if sum(jmp > 0)
                    bnd = bwboundaries(jmp,'noholes');
                    for j = 1:numel(bnd) % for each jump
                        temp = bnd{j}+frst-1; % jump boundaries
                        jb = temp(:,2);
                        mx= max(jb); mn=min(jb);
                        jL = mx-mn+2; % length

                        jSx = tracex(mx+1)-tracex(mn-1); % size
                        jSy = tracey(mx+1)-tracey(mn-1); % size
                        jSz = tracez(mx+1)-tracez(mn-1); % size
                        jsX = jSx/jL;% step
                        jsY = jSy/jL;% step
                        jsZ = jSz/jL;% step
                        jVx = tracex(mn-1)+jsX*(1:jL-1);
                        jVy = tracey(mn-1)+jsY*(1:jL-1);
                        jVz = tracez(mn-1)+jsZ*(1:jL-1);
                        tracex(mn:mx)=jVx;
                        tracey(mn:mx)=jVy;
                        tracez(mn:mx)=jVz;

                    end
                    TraceX(i,:) = tracex;
                    TraceY(i,:) = tracey;
                    TraceZ(i,:) = tracez;
                end                       
                % speed
                xDiff = [tracex(2:end)-tracex(1:end-1) 0];
                yDiff = [tracey(2:end)-tracey(1:end-1) 0];
                zDiff = [tracez(2:end)-tracez(1:end-1) 0];
                traceSpeed = sqrt(xDiff.^2 + yDiff.^2 + zDiff.^2);   

                
                ind = find(tracex~=0);
                frst = min(ind); last = max(ind);
%                 if frst ~= 1, traceSpeed(frst-1) = 0; end
%                 if last ~= Frames ,traceSpeed(last) = 0; end;

                if frst ~= 1
                    traceSpeed(frst-1) = 0; 
                    xDiff(frst-1) = 0; 
                    yDiff(frst-1) = 0; 
                    zDiff(frst-1) = 0; 
                end
                if last ~= Frames 
                    traceSpeed(last) = 0; 
                    xDiff(last) = 0; 
                    yDiff(last) = 0; 
                    zDiff(last) = 0; 
                end;
                xDiff_(i,:) = xDiff;
                yDiff_(i,:) = yDiff;
                zDiff_(i,:) = zDiff;
                TraceSpeed(i,:) = traceSpeed;
                if max(traceSpeed)>15
                    % pause;
                end
            end
            save(traceDataFileNm2,'TraceX','TraceY','TraceZ','TraceINT','TraceSpeed','numCombTraces','cfg')
            save('xyzDiff','xDiff_','yDiff_','zDiff_');
        end
    else
        traceDataFileNm3 = [];
    end
    %% 
    clear mnXdiff mnYdiff
    for i = 1:size(xDiff_,2)
        mnXdiff(i) = mean(xDiff_(xDiff_(:,i)>0,i));
        mnYdiff(i) = mean(yDiff_(yDiff_(:,i)>0,i));
    end
    subplot(1,2,1)
    plot(mean(mnXdiff,1))
    subplot(1,2,2)
    plot(mean(mnYdiff,1))
    %% 

%% PLOT TRACES
    if  1 || ~exist(traceDataFileNm3) && ~isCombDorsalVentral 
        load(traceDataFileNm2);
        disp('preparing plot data...')
        isFewTraces = 0;
        if isFewTraces
            N = 1000;
            TraceX = TraceX(1:N,:);
            TraceY = TraceY(1:N,:);
            TraceZ = TraceZ(1:N,:);
            TraceSpeed = TraceSpeed(1:N,:);
        end
        
        TraceX = full(TraceX);
        TraceY = full(TraceY);
        TraceZ = full(TraceZ);
        TraceSpeed = full(TraceSpeed);
    
        isRemoveSelTraces =0;
        ixDel = [];
        if isRemoveSelTraces % remove corrupt traces
            remDorsal = 4;
            if remDorsal==1 
                Xs = [465];
                Ys = [72];
                dYs = 2;
                dZs =5;
                Zs = [27];
                ixDel = find(...
                    ( (Ys(1)-dYs < TraceY) .* (TraceY<Ys(1)+dYs) ) .* ...
                    ( (Zs(1)-dZs < TraceZ) .* (TraceZ<Zs(1)+dZs ) ) + ...
                    (Xs < TraceX)...
                );
            elseif remDorsal == 0 %remVentral
                Xs = [355];
                Ys = Boy1-[192 374 441 453 530 563]/2; % position of Y stripes
                dYs = 2;
                ixDel = find(...
                    ((Ys(1)-dYs)<TraceY(:)) .* (TraceY(:)<(Ys(1)+dYs)) +...
                    ((Ys(2)-dYs)<TraceY(:)) .* (TraceY(:)<(Ys(2)+dYs)) +...
                    ((Ys(3)-dYs)<TraceY(:)) .* (TraceY(:)<(Ys(3)+dYs)) +...
                    ((Ys(4)-dYs)<TraceY(:)) .* (TraceY(:)<(Ys(4)+dYs)) +...
                    ((Ys(5)-dYs)<TraceY(:)) .* (TraceY(:)<(Ys(5)+dYs)) +...
                    ((Ys(6)-dYs)<TraceY(:)) .* (TraceY(:)<(Ys(6)+dYs)) +...
                    (Xs < TraceX(:))...
                    );
            elseif remDorsal == 3 % showing internal structures
                Xs = [355];
                ixDel = find(...
                    (Xs(1)) > TraceX(:) ...
                    );
            else
                Xs = [300];
                dXs = 60;
                ixDel = find(...
                    ((Xs(1)-dXs) < TraceX(:)) .*  (TraceX(:) < (Xs(1)+dXs) ) ...
                    );
            end
            TraceX(ixDel)=0;
            TraceY(ixDel)=0;
            TraceZ(ixDel)=0;
            TraceSpeed(ixDel)=0;
            ixDel = find(sum(TraceY,2)==0);

            load(traceDataFileNm2);
            isPlotDel = 1;
            if isPlotDel % plots the deleted traces
                TraceX=TraceX(ixDel,:);
                TraceY=TraceY(ixDel,:);
                TraceZ=TraceZ(ixDel,:);
                TraceSpeed=TraceSpeed(ixDel,:);
            else
                TraceX(ixDel,:)=nan;
                TraceY(ixDel,:)=nan;
                TraceZ(ixDel,:)=nan;
                TraceSpeed(ixDel,:)=nan;
            end
        end
        save(traceDataFileNm3,'TraceX','TraceY','TraceZ','TraceINT','TraceSpeed','numCombTraces','ixDel','cfg')
    end
    
    if isCombDorsalVentral
        load DV; % dorsal ventral separation positions
        Fdorsal = dir('dorsal3D\-coeff*');
        Fventral = dir('ventral3D\-coeff*');
        if numel(Fdorsal) > 1 || numel(Fventral) > 1, error('more than one folder'); end;
        fTraceDorsal = dir(['dorsal3D\' Fdorsal.name '\traceData3-coeff*']);
        fTraceVentral = dir(['ventral3D\' Fventral.name '\traceData3-coeff*']);
        
        load(['ventral3D\' Fventral.name '\' fTraceVentral.name]);
        TraceXventral = TraceX;
        TraceYventral = TraceY;
        TraceZventral = TraceZ;
        TraceSpeedventral = TraceSpeed;
        TraceINTventral = TraceINT;
        load(['dorsal3D\' Fdorsal.name '\' fTraceDorsal.name]);
        TraceZ = TraceZ+DV(3)-DV(1)+1; % Dorsal offset 

        TraceX = [TraceX ; TraceXventral];
        TraceY = [TraceY ; TraceYventral];
        TraceZ = [TraceZ ; TraceZventral];
        TraceSpeed = [TraceSpeed ; TraceSpeedventral];
        TraceINT = [TraceINT ; TraceINTventral];
        fPosDataDorsal = dir(['dorsal3D\' Fdorsal.name '\posData-coeff*']);
        fPosDataVentral = dir(['ventral3D\' Fventral.name '\posData-coeff*']);
        load(['dorsal3D\' Fdorsal.name '\' fPosDataDorsal.name],'En1'); %(frames) use the value from tracker function generating trace values
        En1dorsal = En1;
        load(['ventral3D\' Fventral.name '\' fPosDataVentral.name],'En1'); %(frames) use the value from tracker function generating trace values
        En1ventral = En1;
        load(['dorsal3D\' Fdorsal.name '\' fPosDataDorsal.name],'Boy1'); %(frames) use the value from tracker function generating trace values
        imSz = [max(En1dorsal,En1ventral),Boy1];
        load(['dorsal3D\' Fdorsal.name '\' fPosDataDorsal.name],'sptReAppearTime'); %(frames) use the value from tracker function generating trace values
        load(['dorsal3D\' Fdorsal.name '\inputInfo.mat' ],'acqRate_sec','PixelSize','PlaneDist'); %(frames) use the value from tracker function generating trace values
        %save(traceDataFileNm3,'TraceX','TraceY','TraceZ','TraceINT','TraceSpeed')
        
    else
        load(traceDataFileNm3)
    end

    NAN = find(isnan(TraceX));
    TraceX(NAN)=0;
    TraceY(NAN)=0;
    TraceZ(NAN)=0;
    TraceSpeed(NAN)=0;
    
    
    
    imSz = [En1 Boy1];
    
    speedPx2nm = PixelSize/acqRate_sec*1000;
    %pxSzZ = 0.1/6; % um (isotropic)
    pxSzZ = PlaneDist;
    %speedPx2nm=1;
    isHist = 0;
    if isHist
        figure
        hist(TraceSpeed((TraceSpeed(:)>0) )*speedPx2nm,1000);
        title('speed distribution of each trace step')
        xlabel('trace speed (nm/sec)')
        colormap('gray');
    end
    
    % Tracespeed mean
    if 0
        for i = 1:size(TraceSpeed,1) % trace mean
            TraceSpeed(i,TraceSpeed(i,:)>0)= mean(TraceSpeed(i,TraceSpeed(i,:)>0));
        end
    end
    
    if 1
    nAv=3;
        TraceSpeed = full(TraceSpeed);
        for i = 1:size(TraceSpeed,1) % walking average
            TraceSpeed(i,TraceSpeed(i,:)>0)= conv(TraceSpeed(i,TraceSpeed(i,:)>0),ones(nAv,1),'same')/nAv;
        end
    end
    if isHist
        figure
        hist(TraceSpeed((TraceSpeed(:)>0) )*speedPx2nm,1000);
        title('speed distribution of each trace step averaged by neighbours (WAby3)')
        xlabel('trace speed (nm/sec)')
        colormap('gray');
    end
        
    % PLOT
    
    
    isTestColorCoding = 0;
    if isTestColorCoding
        TraceSpeed = 10:64;
        TraceX = TraceSpeed*3;
        TraceY = TraceSpeed*3;
        TraceSpeed(numel(TraceSpeed)+1)=-0.15;
    end
    %TraceZ = max(TraceZ(:))-TraceZ+1;
    % color data
    CplotVecN = size(TraceX,1); % # of traces
    Nframe = size(TraceX,2); % # of frames
    CMinZ = 0; % oclor map in Z
    CMinSpeed = 0; % color in speed
    %speedMax = 3.5;
    if CMinZ 
        Zmin = min(TraceZ(TraceZ~=0));
        Zmax = max(TraceZ(:));
        Zrange = Zmax - Zmin;
        Cplot = round(63*(TraceZ-Zmin) / Zrange)+1;
        CM = colormap('jet');
        ColorMode = 'Z';
    elseif CMinSpeed
        Zmin = min(TraceSpeed(TraceSpeed~=0));
        Zmax = max(TraceSpeed(:));
        Zrange = Zmax - Zmin;
        Cplot = round(63*(TraceSpeed-Zmin) / Zrange)+1;
        CM = colormap('jet');
        ColorMode = 'Speed';
    else % color by trace
        CplotVec = 1: CplotVecN;
        CplotVec = rem(CplotVec,64)+1;
        Cplot = repmat(CplotVec',1,Nframe);
        CM = colormap('lines');
        ColorMode = 'Traces';
    end
    Frames = size(TraceX,2);
    En = Frames;
    [Boy2]=size(TraceX,1);
    
    % select traces
    
    if 0 
        ixSel = 1500:2063;
        %ixSel = 25:50;
        TraceX=TraceX(ixSel,:);
        TraceY=TraceY(ixSel,:);
        TraceZ=TraceZ(ixSel,:);
        TraceSpeed=TraceSpeed(ixSel,:);
        Cplot=Cplot(ixSel,:);
    end
    
    nonZeroIx = find(TraceY>0);
    dispPx = {'*',2}; % plot3k
    minimum = 1;
    maximum = max(Cplot(:));
   
    imgZFout = 'TraceImage';
    magImg = 2;
    mag = 2;
    isScatterPlot = 0; % 1:single image 0: movie
    if ~isScatterPlot, mag = 1; end
    imgZFout = 'TraceImage';
    TraceY = imSz(2)-TraceY+1;
    TraceX = TraceX*mag;
    TraceY = TraceY*mag;
    
    figSz(1) = imSz(1)*magImg;
    figSz(2) = imSz(2)*magImg;
    colormap('gray');

    zMax = max(TraceZ(:));
    if zMax == 0, zMax = 1; TraceZ = TraceZ+1; end;

    % find the frames where the traces disappear
    hat = ones(1,sptReAppearTime); 
    dspTrcFrm = zeros(size(TraceX,1),1);
    for i = 1:size(TraceX,1) % all traces
        trace = full(TraceX(i,:));
        trFrm = conv(trace,hat,'same'); % active frames of the trace
        trFrm(1:find(trFrm>0,1)) = 1;
        lastFrm = find(trFrm==0,1); 
        if ~isempty(lastFrm)
            dspTrcFrm(i) = lastFrm;
        end
    end

    imgZFout = 'TraceImageStack';
    imgZFoutScat = 'TraceImage0';
    fOut = sprintf('%s_%s.tif',imgZFout,ColorMode);
    fOutScat = sprintf('%s_%s.tif',imgZFoutScat,ColorMode);
    hQ = 0; %hImg = image; 
    lastX = TraceX(:,1); 
    lastY = TraceY(:,1);
    lastZ = TraceZ(:,1);
    tit = 'image';
    m = imSz(1)*mag; n = imSz(2)*mag;
    pos=get(0,'ScreenSize');
    pos=pos(3:4) - [m n-35];
    %fig=figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos/2 m n]);
    fig=figure('DoubleBuffer','on','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos/2 m n]);
    axe=axes('Parent',fig,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren','XLim',0.5+[0 m],'YLim',0.5+[0 n]);

    %% TIF description
    if isCombDorsalVentral
        tiffDescription = 'dorsal ventral combined';
    else
        tiffDescriptionImg = sprintf('\nImage Info: \nPixelSize:%.03fum \nframeTime:%.03f seconds  ',...
            PixelSize/magImg,frameTime);
        tiffDescriptionData = sprintf('\nData Info: \nPixelSize:%.02fum \nPlaneDist:%.02fum \nStackNum:%i \nFrames:%i ',...
            PixelSize,PlaneDist,StackNum,Frames);
        tiffDescriptionDetection = sprintf('\nDetection: \ngausWin:%.02f pixels \ngausSg:%.02f pixels \nlap (3by3):%s \nCoeff:%.02f ',...
            gausWin,gausSg,sprintf('%i ',lap),Coeff);
        tiffDescriptionLocalization = sprintf('\nLocalization: \nWindowSize:%.02f ',...
            WindowSize);
        tiffDescriptionTrace = sprintf('\nTracing: \nsptJmp:%.02f pixels \nsptJmpCombination:%.02f pixels \nsptReAppearTime:%.02f frames \nminTraceLength:%.02f frames \nNumber of Trace Combination:%i',...
            sptJmp,sptJmpCombination,sptReAppearTime,minTraceLength,numCombTraces);
        tiffDescription = sprintf('\n %s\n %s\n %s\n %s\n %s',...
            tiffDescriptionImg,tiffDescriptionData,tiffDescriptionDetection,tiffDescriptionLocalization,tiffDescriptionTrace);
    end
    %% PLOTTING    colormap('jet')
    if isScatterPlot
        xx = TraceX(TraceSpeed>0);
        yy = TraceY(TraceSpeed>0);
        zz = TraceZ(TraceSpeed>0);
        cc = CM(Cplot(TraceSpeed>0),:);
        isShowTraceLine = 1; % Lines or Scatter Plot
        if isShowTraceLine
            TraceX(TraceX==0)=nan;
            TraceY(TraceX==0)=nan;
            line(TraceX',TraceY')
        else
            scatter3(xx,yy,zz,4,cc,'*');
        end        
        view(0,90);
        gry = 0.5;
        set(gcf,'Color',[1 1 1]*gry);
        imgFig = getframe(gcf);
        imgOut = imgFig.cdata;
        figPos = get(gcf,'Position');
        %imgOut = imgOut((figPos(4)-m)/2+1:(figPos(4)+m)/2,(figPos(3)-n)/2+1:(figPos(3)+n)/2,:);
        imgOut = imgOut((figPos(4)-n)/2+1:(figPos(4)+n)/2,(figPos(3)-m)/2+1:(figPos(3)+m)/2,:);
        imwrite(imgOut,fOutScat,'Compression', 'none') 
        colorbarFout = sprintf('%s_colorbar.tif',fOutScat(1:end-4));
        setTiffDescription(fOutScat,tiffDescription)   
    else
        colormap('gray')
        
        img2D = imread(fnameOrigin,1); 
        img2D = flipud(img2D(yy1:yy2,xx1:xx2));
        hImg = imagesc(img2D,'Parent',axe); %axis image; 
        set(gca,'position',[0 0 1 1]); set(gcf,'position',[0 0 figSz(1) figSz(2)]);
        
        imgFig = getframe(gcf);
        imgOut = imgFig.cdata;
        figPos = get(gcf,'Position');
        if exist([imgZFout '.tif'])
            delete([imgZFout '.tif']);
        end
        imwrite(imgOut,fOut,'Compression', 'none') 

        for ixFrm = 1:En-1
            img2D = imread(fnameOrigin,ixFrm+1); 
            img2D = flipud(img2D(yy1:yy2,xx1:xx2));
            

            hImg = imagesc(img2D,'Parent',axe); %axis image; 
            %set(gca,'position',[0 0 1 1]); set(gcf,'position',[0 0 figSz(1) figSz(2)]);
            ixTrc = find(TraceX(:,ixFrm+1)>0);
            currX = TraceX(:,ixFrm+1);
            currY = TraceY(:,ixFrm+1);
            currZ = TraceZ(:,ixFrm+1);
            uistack(hImg,'bottom'); hold on
            ixTrc2 = find(lastX(ixTrc)>0);
            ix = ixTrc(ixTrc2); % indices of traces tracked in the current frame    
            dspTrcIx = find(ixFrm+1 == dspTrcFrm); % index for dissappearing traces
            showTrace =1;
            if showTrace
                hQdel = hQ(dspTrcIx,:); % handles for discont. traces
                delete(hQdel(find(hQdel~=0)));    % remove the traces of the discontinued traces
            end

            for i = 1:round(length(ix)) % for each trace 
                iL = ix(i);  % index for each line
                hQ(iL,ixFrm)=quiver(lastX(iL),lastY(iL),currX(iL)-lastX(iL),currY(iL)-lastY(iL),'Color',CM(Cplot(iL,ixFrm),:));
                if ~showTrace
                    adjust_quiver_arrowhead_size(hQ(iL,ixFrm-1),5)
                end
            end
            % print images
            imgFig = getframe(gcf);
            imgOut = imgFig.cdata;
            figPos = get(gcf,'Position');
            %imgOut = imgOut((figPos(4)-m)/2+1:(figPos(4)+m)/2,(figPos(3)-n)/2+1:(figPos(3)+n)/2,:);
            %imgOut = imgOut((figPos(4)-n)/2+1:(figPos(4)+n)/2,(figPos(3)-m)/2+1:(figPos(3)+m)/2,:);

            imwrite(imgOut,fOut,'WriteMode','append','Compression', 'none') 

            delete(hImg)
            %export_fig(imgZFout,'-append');
            fprintf('frame %i/%i \n',ixFrm+1,En);
            lastX = currX; 
            lastY = currY;
            lastZ = currZ;
        end
        setTiffDescription(fOut,tiffDescription)        
    end
    return
    %% colorbar
    if CMinZ 
        plotColorBar(TraceZ*pxSzZ,7,'Z [\mum]',0.5); % 
    elseif CMinSpeed
        plotColorBar(TraceSpeed*speedPx2nm,5,'speed [nm/sec]',0.5);
    end
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'Color',[1 1 1]);
    imgFig = getframe(gcf);
    imgOut = imgFig.cdata;
    imwrite(imgOut,colorbarFout,'Compression', 'none') 
    
%    [TraceX(1135,85:90)' TraceY(1135,85:90)' TraceZ(1135,85:90)']


return
    if isWrite2XLS
        load(traceDataFileNm)
        % remove nan rows
        ixnan = find(isnan(TraceX(:,1)));
        for i = 1:numel(ixnan)
            TraceX = TraceX([1:ixnan(i)-1 ixnan(i)+1:end],:);
            TraceY = TraceY([1:ixnan(i)-1 ixnan(i)+1:end],:);
            TraceZ = TraceZ([1:ixnan(i)-1 ixnan(i)+1:end],:);
            ixnan = ixnan - 1;
        end
        % write to excel
        fname = 'TraceX.xls';
        xlswrite(fname,TraceX)
        fname = 'TraceY.xls';
        xlswrite(fname,TraceY)
        fname = 'TraceZ.xls';
        xlswrite(fname,TraceZ)
    end
    
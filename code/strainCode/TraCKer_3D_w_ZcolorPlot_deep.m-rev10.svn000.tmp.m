        clear; close all;
    F = findall(0,'type','figure','tag','TMWWaitbar'); delete(F)
    isFindRecruitment = 0;
    PlaneDist = 3;
    PixelSize = 1;
    isCropXY = 0;
    isDebugFilt = 0;
    isWalkAvFilt = 0;
    isWAby3 = 0; %  walking averga by 3 
    frstFrm = 1000;
    frstFrm = 1;
    
    % the min. distance between spots of a trace in consecutive frames (spot jump)
    sptJmp = 6; %(pixels)
    sptReAppearTime = 2; %(frames)
    minTraceLength = 4;
    minTraceLength = 2;

    isCombTraces = 0;
    isGausFit = 0; % Gaussian localization
    isGausFitLocalize = 1; % gaussian fit filtering and localization
    StackNum = 1;
    WindowSize = 5; 
    BigWindowSize=WindowSize+4;
    %BigWindowSize=WindowSize;
    
    % selects the background for laser intensity normalization
    fnLog = 'logFile.txt';
    fidLog = fopen(fnLog,'a+');

    %READ the ORIGINAL file
    if exist('fname.mat')
        load fname
        if ~exist(fname)
            fname = sprintf('..\\%s',fname);
        end
    end
    fTif = dir('*.tif');
    for i = 1:length(fTif)
        if ~isempty(strfind(fTif(i).name,'_000.tif')) && isempty(strfind(fTif(i).name,'MAX'))
            fnameBck = fTif(i).name;
        end
    end
    imgFrst = imread(fname);
    [Boy1,En1]=size(imgFrst);
    if isCropXY
        xx1 = 300; xx2 = 700; % crop
        yy1 = 300; yy2 = 700; % crop
    else
        xx1 = 1; xx2 = En1; % crop
        yy1 = 1; yy2 = Boy1; % crop
    end
    [Boy1,En1]=size(imgFrst(yy1:yy2,xx1:xx2));
    sizeImg=[Boy1,En1];
    %imageInfo=imfinfo([fname(1:end-4) '_' fname(end-3:end)]);
    imageInfo=imfinfo(fname);
    Frames=length(imageInfo);
    Frames=5000;
    
    imageInfo = imageInfo(1);
    %Frames = 20; % # of frames to be analyzed

    coeffMat = dir('coeff*.mat');
    Coeff = 0;
    
    coeffFound = 0;
    if ~isempty(coeffMat)
        if numel(coeffMat) > 1
            for i = 1:numel(coeffMat)
                if strcmp(fname(end-6:end-4),coeffMat(i).name(end-6:end-4))
                    coeffFound = coeffFound+1;
                    iCoeff = i;
                end
                if coeffFound>1
                    error('multiple coeff files with the same number')
                end
            end
            if coeffFound == 1
                coeffMat = coeffMat(iCoeff).name;
            end
        elseif length(coeffMat.name) == 13 && strcmp(fname(end-6:end-4),coeffMat.name(end-6:end-4))
            coeffMat = coeffMat.name;
            coeffFound = 1;
        else
            display('No Coeff file was found associated with the movie file');
        end
        if coeffFound == 1
            load(coeffMat)
            display('loading coefficient')
        end
    else
        display('No Coeff file was found');
    end
    assignFileNames

    %% FIND X & Y
    if ~exist(spotWinFileNm) 
        if ~exist('BWbckgrnd.mat'), selBackgrnd; else, load('BWbckgrnd.mat'), end;
        topHat4by4 = [-0.0833   -0.0833   -0.0833   -0.0833;-0.0833    0.2500    0.2500   -0.0833;-0.0833    0.2500    0.2500   -0.0833;-0.0833   -0.0833   -0.0833   -0.0833]; % topHat
        topHat5by5 = [-0.0625   -0.0625   -0.0625   -0.0625   -0.0625;-0.0625    0.1111    0.1111    0.1111   -0.0625;-0.0625    0.1111    0.1111    0.1111   -0.0625;-0.0625    0.1111    0.1111    0.1111   -0.0625;-0.0625   -0.0625   -0.0625   -0.0625   -0.0625];
        topHat5by5 = [-0.0649   -0.0649   -0.0649   -0.0649   -0.0649; -0.0649    0.1087    0.1087    0.1087   -0.0649;-0.0649    0.1087    0.1087    0.1087   -0.0649;-0.0649    0.1087    0.1087    0.1087   -0.0649;-0.0649   -0.0649   -0.0649   -0.0649   -0.0649];
        topHat7by7 = [-0.0417   -0.0417   -0.0417   -0.0417   -0.0417   -0.0417   -0.0417;-0.0417    0.0400    0.0400    0.0400    0.0400    0.0400   -0.0417;-0.0417    0.0400    0.0400    0.0400    0.0400    0.0400   -0.0417;-0.0417    0.0400    0.0400    0.0400    0.0400    0.0400   -0.0417;-0.0417    0.0400    0.0400    0.0400    0.0400    0.0400   -0.0417;-0.0417    0.0400    0.0400    0.0400    0.0400    0.0400   -0.0417;-0.0417   -0.0417   -0.0417   -0.0417   -0.0417   -0.0417   -0.0417];
        lap=[-1,-1,-1;-1,8,-1;-1,-1,-1];
        gaus1 = ones(2)/4;
        gaus2 = fspecial('gaus',3,1);
        gaus2 = gaus2/sum(gaus2(:));
        backgrnd = ones(9);
        gausHat1 = -backgrnd/(numel(backgrnd)-numel(gaus1));
        gausHat2 = -backgrnd/(numel(backgrnd)-numel(gaus2));
        gausHat1(4:5,4:5)=gaus1;
        gausHat2(4:6,4:6)=gaus2;
        elev=1/(81-25);
        gausHat1=fspecial('gaussian', 5, 1);

        gausHat1 = gausHat1 + elev;
        gausHatPad1 = padarray(gausHat1,[2 2]);
        gausHatPad1 = gausHatPad1 - elev;
        gausHat1 = gausHatPad1;
        gausHat2 = gausHat1;
        
        if ~coeffFound
            % generate preview image
            nTile = 0;
            genPreviewImg;  % generate the image SHOW
            
            imSize = size(SHOW);
            m = imSize(1); n = imSize(2);
            tit = 'determine the coefficient for background filtering';
            hFigDisp=figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256));
            set(gcf,'units','normalized','outerposition',[0 0 1 1])
            kk =0.1;
            axe=axes('Parent',hFigDisp,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','outerposition',[0 0 1 1],'Nextplot','replacechildren','XLim',0.5+[0 m],'YLim',0.5+[0 n]);

            %Scale=[1:1:10000];
            mx = max([PREdataFiltTiled1(:); PREdataFiltTiled2(:)]);
            Scale=[1:1:mx];
            tit2 = 'controls';
            hFigCont = figure('DoubleBuffer','on','Menubar','none','Name',tit2,'NumberTitle','off','Colormap',gray(256));
            figPos = [1921 -461 1080 1844];
            set(hFigCont,'Position', figPos);

            %draw histogram
            axe=axes('Parent',hFigCont ,'Visible','off');
            hist([0 ;double(SHOW(SHOW>0))],1000); histSHOW = hist([0; double(SHOW(SHOW>0))],1000);
            [ maxHistDisp histPeakPos] = max(histSHOW); ytick = get(gca,'YTick'); maxHistDisp = ytick(end);
            
            bt = 40;
            hTextCoeff = uicontrol('style','text','String','coefficient','Position',[20 bt+20 160 15]);
            hCoeff = uicontrol('style','slider','units','pixel','position',[20 bt 300 20]); 
            hTextImgMin = uicontrol('style','text','String','min','Position',[20 bt+60 160 15]);
            hImgMin = uicontrol('style','slider','units','pixel','position',[20 bt+40 300 20]); 
            hTextImgMax = uicontrol('style','text','String','max','Position',[20 bt+100 160 15]);
            hImgMax = uicontrol('style','slider','units','pixel','position',[20 bt+80 300 20]); 
            hTextImgGamma = uicontrol('style','text','String','gamma','Position',[20 bt+140 160 15]);
            hImgGamma = uicontrol('style','slider','units','pixel','position',[20 bt+120 300 20]); 
            imgMax_=double(max(SHOW(:)));
            filtMax = double(max(max([PREdataFiltTiled1;PREdataFiltTiled2])));
            Coeff = filtMax/5;
            q=0; 
            CoeffSlider = double(double(Coeff)/filtMax);
            imgMinSlider = 0;  %double(histPeakPos)/1000
            imgMin = imgMax_*imgMinSlider;
            imgMaxSlider = 1; imgMax = imgMax_*imgMaxSlider;
            imgGammaSlider = 0.1; imgGamma = 1.5; %  minG = 0.2; maxG = 8.2;
            set(hCoeff,'Value',CoeffSlider); 
            set(hImgMin,'Value',imgMinSlider); 
            set(hImgMax,'Value',imgMaxSlider); 
            set(hImgGamma,'Value',imgGammaSlider);  
            
            hLineMin = line([imgMin imgMin],[0 maxHistDisp ]); set(hLineMin,'LineWidth',2); set(hLineMin,'Color',[1 0 1])
            hLineMax = line([imgMax imgMax],[0 maxHistDisp ]); set(hLineMax,'LineWidth',2); set(hLineMin,'Color',[0 1 1])
            nP = 1000*(imgMaxSlider - imgMinSlider); 
            gammaLine = 1:nP;
            gamma = gammaLine.^imgGamma;
            gamma = gamma/max(gamma)*maxHistDisp ;
            imgRange = imgMax-imgMin;
            
            hold; hLineGamma = plot(imgMin+gammaLine/1000*filtMax,gamma); hold;
            BINARprev =0;BINARprev2=0;
            while q == 0 
                if exist('listen1'), delete(listen1);delete(listen2);delete(listen3);delete(listen4); end;
                
                CoeffTiled = BckTiled*Coeff/bck(1);
                %colormap(gray);
                % find peaks
                PREdataFiltTiledDiv=PREdataFiltTiled1./CoeffTiled;
                BINAR = im2bw(PREdataFiltTiledDiv,1);
                DINAR1=uint16(BINAR).*PREdataTiled; %DINAR1 = 0;
                PREdataFiltTiledDiv=PREdataFiltTiled2./CoeffTiled;
                BINAR = im2bw(PREdataFiltTiledDiv,1);
                DINAR2=uint16(BINAR).*PREdataTiled; 
                BWSHOW=imregionalmax(DINAR1+DINAR2, 8);
                [y,x,v]=find(BWSHOW==1);
                %save('filt1','temp1','BWSHOW','BINAR','PREdataFiltTiledDiv','PREdataFiltTiled1','DINAR1','Coeff')
    
                figure(hFigDisp);
                set(gcf,'Name',sprintf('frames:%i - %i',1+nTile*nFrTile,nFrTile+nTile*nFrTile ))
                %plot
                img=SHOW;scaleImg;  % scale image
                if isDebugFilt, subplot(1,4,1); end
                hImg = imagesc(img); axis image;
                hold on;
                hscat = scatter(x,y,22,'o');hold off;
                axis image;

                if isDebugFilt
                    subplot(1,4,2);
                    imagesc(PREdataFiltTiled1); axis image;
                    hold on;
                    hscat2 = scatter(x,y,22,'o');hold off;
                    subplot(1,4,3);
                    imagesc(PREdataFiltTiled2); axis image;
                    hold on;
                    hscat3 = scatter(x,y,22,'o');hold off;
                    subplot(1,4,4);
                    imagesc(BWSHOW); axis image; 
                else
                    hscat2 = [];hscat3 = [];
                end
                
                BINARprev = 0; BINARprev2 = 0;
                figure(hFigCont)
                listen1 = addlistener(hCoeff,'ActionEvent',@(hObject, event) updCoeff(hObject, event,hFigDisp,hTextCoeff,hscat,hscat2,hscat3,bck,sizeImg,nRow,PREdataTiled,PREdataFiltTiled1,PREdataFiltTiled2,filtMax,isDebugFilt)); 
                listen2 = addlistener(hImgMin,'ActionEvent',@(hObject, event) updImgMin(hObject, event,nTileXY,hTextImgMin,hImg,SHOW,imgMax_, hImgMax, hImgGamma,hLineMin,hLineGamma));
                listen3 = addlistener(hImgMax,'ActionEvent',@(hObject, event) updImgMax(hObject, event,nTileXY, hTextImgMax,hImg,SHOW,imgMax_,hImgMin, hImgGamma,hLineMax,hLineGamma));
                listen4 = addlistener(hImgGamma,'ActionEvent',@(hObject, event) updImgGamma(hObject, event,nTileXY, hTextImgGamma,hImg,SHOW,imgMax_,hImgMin, hImgMax,hLineGamma));
                
                figure(hFigDisp);
                btn = 0;
                
                while btn == 0
                    btn = waitforbuttonpress;
                    k = get(hFigDisp,'CurrentCharacter');
                end
                CoeffSlider = get(hCoeff,'Value'); Coeff = CoeffSlider*filtMax;
                upd = 0; isNewTileImage = 0;
                switch lower(k)
                    case 's'
                        CoeffSlider = CoeffSlider - 0.001; upd = 1;
                    case 'd'
                        CoeffSlider = CoeffSlider + 0.001; upd = 1;
                    case 'w'
                        CoeffSlider = CoeffSlider - 0.01; upd = 1;
                    case 'e'
                        CoeffSlider = CoeffSlider + 0.01; upd = 1;
                    case 'q'
                        q = 1;
                    case 'z'
                        return;
                    case 'r'
                        nTile = nTile - 1; isNewTileImage = 1;
                    case 't'
                        nTile = nTile + 1;  isNewTileImage = 1;
                    case 'f'
                        nTile = nTile - 10;  isNewTileImage = 1;
                    case 'g' 
                        nTile = nTile + 10; isNewTileImage = 1;
                end
                if isNewTileImage
                    genPreviewImg; % update image
                end
                if upd == 1
                    set(hCoeff,'Value',CoeffSlider); 
                    Coeff = CoeffSlider*filtMax;
                end
                set(hTextCoeff,'String',sprintf('coefficient: %i',round(Coeff)));
            end
            save_Coeff;
            return
            
        end
        close;
        assignFileNames
        
        fnCrop=dir('bckgrnd_*.txt');  % e.g. bckgrnd_12X76Y51x21.txt
        if ~isempty(fnCrop)
            fnameCrop = fnCrop.name;
            nmCrop = fnameCrop(9:end-4);
            ixX = find(nmCrop=='X');
            ixY = find(nmCrop=='Y');
            ixx = find(nmCrop=='x');
            xCr = str2num(nmCrop(1:ixX-1));
            yCr = str2num(nmCrop(ixX+1:ixY-1));
            szXcr = str2num(nmCrop(ixY+1:ixx-1));
            szYcr = str2num(nmCrop(ixx+1:end));
        end
        
        
        if isWalkAvFilt
            for j=1:Frames % read images
                temp = imread(fname,j);
                J(:,:,j) = temp(yy1:yy2,xx1:xx2);
                %bckgrndImg(:,:,j) = J(yCr+1:yCr+szYcr,xCr+1:xCr+szXcr,j);
                %bckgrnd(j) = sum(sum(J(yCr+1:yCr+szYcr,xCr+1:xCr+szXcr,j)));
            end
            J=uint16(J);
        end
        
        isLaserFluct = 0;
        if isLaserFluct
            [szY,szX] = size(bckgrndImg);
            aaa= reshape(bckgrndImg,[szX*szY,Frames]);
            medImage = median(aaa,1)
            bckgrndMx = max(bckgrnd);
            bckgrndNorm = bckgrnd/bckgrndMx;
            
        end
        
        
        if isWalkAvFilt
            for j=1:Frames-2 % filter images
                if isWAby3 % walking averga by 3 isWAby3
                    imgWA = (double(J(:,:,j))+double(J(:,:,j+1))+double(J(:,:,j+2)))/3; 
                else
                    imgWA = (double(J(:,:,j))+double(J(:,:,j+1)))/2; 
                end
                img = imgWA;
            end
        end
        
        % read Z data
        zFile = 'Z.txt';
        if exist(zFile)
            fidZ = fopen(zFile);
            tLine=fgetl(fidZ);
            lineNum=1;
            while ischar(tLine)
                Zread(lineNum)=double(str2num(tLine));
                tLine=fgetl(fidZ);
                lineNum=lineNum+1;
            end
            Zread = Zread-min(Zread);
        end
        

        h = waitbar(0,'3D localization...');
        %% find 3D intensity
        sp = 1;
        nSp = 0;
        dinPrev2 = 0;dinPrev = 0;
        
        for k=1:Frames-frstFrm+1
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
            temp = imread(fname,k+frstFrm-1); 
            IMG = double(temp(yy1:yy2,xx1:xx2));
            
            temp2 = imread(fnameBck,k+frstFrm-1); 
            temp2 = temp2.*uint16(BWbckgrnd);   
            bck(k) = median(median(double(temp2(temp2>0))));
            
            % background
            
            
            if k == 1, bckFirst=bck(1); end;
            
            CoeffNorm(k) = bck(k)*Coeff/bckFirst;
            CoeffNorm(k) = round(CoeffNorm(k));
            
            dataFilt = imfilter(IMG,gausHat1,'symmetric');
            dataFilt = imfilter(dataFilt,lap,'symmetric');
            %dataFilt2(:,:,j) = imfilter(IMG,gausHat2,'symmetric');
            %dataFilt2(:,:,j) = imfilter(dataFilt2(:,:,j),lap,'symmetric');
            
            nxt2ndFrm = k+2; if k+2>Frames,nxt2ndFrm=Frames;end;
            % find peaks
%            CoeffNorm = round(Coeff*bckgrndNorm);
             
            %dataFilt=uint16(dataFilt);
            dataFiltDiv = dataFilt/CoeffNorm(k);
            bin = im2bw(dataFiltDiv,1);
            din1 = uint16(bin).*uint16(IMG);
%             dataFiltDiv = dataFilt2/CoeffNorm;
%             bin = im2bw(dataFiltDiv,1);
%             din2 = uint16(bin).*IMG; 
din2 = 0;
            if isempty(find( din1 + din2 > 0)), continue; end;
            BW = imregionalmax(din1 + din2, 8);
            %save('filt2','IMG','BW','bin','dataFiltDiv','dataFilt','din1')
            [B,L] = bwboundaries(BW,'noholes');

            q=0; 
            %DEFINE Size
            %[En,Boy]=size(IMG);
            nSp = nSp + length(B);
            for m=1:length(B) % for each spot
                c=cell2mat(B(m));
                %csize=(max(c(:,1))-min(c(:,1)))*(max(c(:,2))-min(c(:,2)));
                q=q+1;
                PyBW=uint16(mean(c(:,1)));
                PxBW=uint16(mean(c(:,2)));
                
                % ignore the spots on the edge
                if ((PxBW-(BigWindowSize+1)/2)<1 || (PyBW-(BigWindowSize+1)/2)<1 || (PxBW+(BigWindowSize+1)/2)>En1 || (PyBW+(BigWindowSize+1)/2)>Boy1)
                    q = q-1;
                    continue
                end
                %DEFINE Big Window
                yy = PyBW-(BigWindowSize+1)/2;
                y1 = yy + 1; y2 = yy + BigWindowSize;
                xx = PxBW-(BigWindowSize+1)/2;
                x1 = xx + 1; x2 = xx + BigWindowSize;
                BigWindow=IMG(y1:y2,x1:x2);   
                
                spotWin(:,:,sp) = BigWindow*bckFirst/bck(k);  
                xBW(q,k)=PxBW; 
                yBW(q,k)=PyBW;
                if exist('Zread')
                    zBW(q,k)=Zread(k);
                else
                    zBW=zeros(size(xBW));
                end
                sp = sp + 1;
            end
            waitbar(k/Frames)
            if ~rem((k-1),100) 
                time100frame = 0;
                if k > 1
                    time100frame = toc; 
                    %disp(sprintf('time for 100 frames : %.02f\n',time100frame));
                end;
                tic;
                
            end
        end
        close(h);
        save(spotWinFileNm,'spotWin','xBW','yBW','zBW'); 
        clear X Y Z INT spotWin NBINs;
        clear JF J bins JF1 JF2 JF0 yBW xBW Yx Xc BACK
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
            copyfile(coeffMat,strcat(folderNm,'\',nmCoeff));
            %copyfile(fname,strcat(folderNm,'\',fname));
            if exist('spotSel.mat'), copyfile('spotSel.mat',strcat(folderNm,'\spotSel.mat')); end;
            if exist('spotSelVal.mat'), copyfile('spotSelVal.mat',strcat(folderNm,'\spotSelVal.mat')); end
            
            movefile(posDataFileNm,strcat(folderNm,'\',posDataFileNm));
            movefile(spotWinFileNm,strcat(folderNm,'\',spotWinFileNm));  
            %copyfile('fname.mat',strcat(folderNm,'\fname.mat'));          
            cd(folderNm)
            fname = ['..\' fname ];save('fname','fname');
        end
        delete(xyzDataGausFileNm)
    end
    %return
    
    %% DOUBLE GAUS FIT
    if isGausFitLocalize && ~exist(xyzDataGausFileNm)
        % load data
        fnSpot = rdir('spotWin*');
        load(char(fnSpot.name))
        spotWin = double(spotWin);
        % load parameters

        % parameters
        load spotSel;
        %intP1G = 3; intP = 3; intP2G = 3; % needs to be odd
        %tP2G_1 = 1; % tolerance position
        %tP2G_2 = 0.2;        
        %save_spotSel
        % sg shoud be around 220 nm
        % pixel size : 106nm
        sg = FWHM/(2*sqrt(2*log(2)))*sgBlur; % gauss sigma/(2*sqrt(2*log(2)))*sgBlur; % gauss sigma
        sg= 2.13/1.6;
        sg2 = 5/2+1; sr2 = inf;
        minmaxIntTol = 0.6;
        tol = minmaxIntTol + 1;
        sr = sr0+1;
        sg1G = sg*sigma1Gvs2G;
        sr1G = sr0*sigma1Gvs2G+1;
        sg2by2 = 0.77;

        % define output XY array
        [iyTrack ixTrack] = find(xBW>0); % positions of spot data in the X array
        Xgaus = zeros(size(xBW));
        Ygaus = zeros(size(xBW));
        Xgaus2by2 = zeros(size(xBW));
        Ygaus2by2 = zeros(size(xBW));

        % crop 
        x3x1=4;y3y1=x3x1;x3x2=6;y3y2=x3x2;
        x5x1=3;y5y1=x5x1;x5x2=7;y5y2=x5x2;
        
%        x3x1=2;y3y1=x3x1;x3x2=4;y3y2=x3x2;
%        x5x1=1;y5y1=x5x1;x5x2=5;y5y2=x5x2;
        Window5by5 = spotWin(y5y1:y5y2,x5x1:x5x2,1);

        %define fit parameters
        [n,m]=size(Window5by5);%assumes that I is a nxm matrix
        px0 = double(ceil(m/2));
        
        c01G5by5 = double([0 1 px0 sg1G px0 sg1G]); %start-guess here
        fun1 = @(c,x) c(1)+c(2)*exp(-intP1G*((x(:,1)-c(3))/c(4)/sqrt(2)).^2-intP1G*((x(:,2)-c(5))/c(6)/sqrt(2)).^2);
        lb1G5by5 = double([0 1/tol 1 sg1G/sr1G 1 sg1G/sr1G]); 
        ub1G5by5 = double([inf tol m sg1G*sr1G m sg1G*sr1G]);
        % 2 gaussian
        TBA = 0;
        c02G5by5([4,6,9,11]) = double([sg sg sg2 sg2]); %start-guess here
        lb2G5by5([4,6,9,11]) = double([sg/sr sg/sr sg2 sg2]);  % [back peak1 x1 sigX1 y1 sigY1 peak2 x2 sigX2 y2 sigY2]
        ub2G5by5([4,6,9,11]) = double([sg*sr sg*sr sg2*sr2 sg2*sr2]); 
        fun2 = @(c,x) c(1)+c(2)*exp(-intP*((x(:,1)-c(3))/c(4)/sqrt(2)).^2-intP*((x(:,2)-c(5))/c(6)/sqrt(2)).^2) + c(7)*exp(-intP*((x(:,1)-c(8))/c(9)/sqrt(2)).^2-intP*((x(:,2)-c(10))/c(11)/sqrt(2)).^2);
        fun2_2 = @(c,x) c(1)+c(2)*exp(-intP2G*((x(:,1)-c(3))/c(4)/sqrt(2)).^2-intP2G*((x(:,2)-c(5))/c(6)/sqrt(2)).^2) + c(7)*exp(-intP2G*((x(:,1)-c(8))/c(9)/sqrt(2)).^2-intP2G*((x(:,2)-c(10))/c(11)/sqrt(2)).^2);

        % gaus2by2 fit grid
        X_ = 0.5:0.2:2.5;
        [XX,YY]=meshgrid(X_,X_);
        XY = [XX(:) YY(:)];
        
        % gaus fit grid
        [XX,YY]=meshgrid(1:n,1:m);%your x-y coordinates
        x(:,1)=XX(:); % x= first column
        x(:,2)=YY(:); % y= second column
        
        % gaus display grid
        [n2,m2]=size(spotWin(:,:,1)); %assumes that I is a nxm matrix
        [XX,YY]=meshgrid(1:n2,1:m2); %your x-y coordinates
        xSHOW(:,1)=XX(:); % x= first column
        xSHOW(:,2)=YY(:); % y= second column
        
        h = waitbar(0,'Gaus fit localization ....');
        nSpots = size(spotWin,3);
        %nSpots = 200;
        nIntPeakAtEdge = 1;
        frstSpot =1;
        sp = frstSpot;
        errCenter = cell(1,nSpots);
        isFit5by5 = cell(1,nSpots); 
        for i = frstSpot:nSpots % every loop

            is1GinCenter(sp) = 1;
            isMaxinCenter(sp) = 1;
            isSpot(sp) = 1;
            %% gauss fit 1/2
            % min & max values
            Window5by5 = double(spotWin(y5y1:y5y2,x5x1:x5x2,i));
            [py5Ymax,px5Xmax] = find(Window5by5==max(Window5by5(:)));
            if numel(py5Ymax) > 1
                R_ = (py5Ymax-3).^2+(px5Xmax-3).^2;
                [vv ix] = sort(R_);
                py5Ymax=py5Ymax(ix(1)); px5Xmax=px5Xmax(ix(1));
            end
            intMax(sp) = max(max(Window5by5)); % peak intensity
            Window5by5 = Window5by5/intMax(sp);
            backGrnd=median(median(Window5by5)); % background
            intPeak = 1 - backGrnd;
            backGrnd_(sp)=backGrnd;
            Window5by5 = Window5by5/intPeak; % spot peak is normalized to one
            backGrnd = backGrnd/intPeak;
            Idat(:,:,sp) = Window5by5;
            if py5Ymax < 2 || py5Ymax > 4 || px5Xmax < 2 || px5Xmax > 4 % max. peak is not in the 3 by 3 window
                errCenter(sp) = {'Max peak is not in the center'};
                intPeakAtEdge(sp) = intPeak;
                nIntPeakAtEdge = nIntPeakAtEdge + 1;
                intMax(sp) = max(max(double(spotWin(y3y1:y3y2,x3x1:x3x2,i)))) - backGrnd;
                isMaxinCenter(sp) = 0;
            end



            intPeak_(sp) = intPeak; 
            % peak pos.
            [cc1(:,sp),err1G5by5(sp),err1res_,EXITFLAG(1)] = gausFit(Window5by5,fun1,c01G5by5,lb1G5by5,ub1G5by5);

            %% err2by2
            pxX = cc1(3,sp);
            pxY = cc1(5,sp);
            pxXc = round(cc1(3,sp));
            pxYc = round(cc1(5,sp));
            meanSpotPeak(sp)=nan;
            isTryFit(sp) = 1; 
            PX2by2(sp) = nan; PY2by2(sp) = nan;
            if pxXc==1 || pxYc==1 || pxXc==m || pxYc==m  % 1Gfit peak is on the edge
                is1GinCenter(sp) = 0;
                isFit5by5(sp) = {'1Gfit peak not in the center'};
                isTryFit(sp) = 0;
            else % find 2by2 and 3by3 info
                data3by3 = Idat(pxYc-1:pxYc+1,pxXc-1:pxXc+1,sp); 
                % find the 4 peak pixels
                iserr2by2_ = [];
                dat=data3by3;
                clear ixX ixY;
                datCenter = dat(2,2); dat(2,2)=0;
                [v iX]=max(dat(:)); dat(iX) = 0;
                ixY(1) = 2; ixX(1) = 2;
                [ixY(2) ixX(2)]=ind2sub(size(dat),iX);
                [v iX]=max(dat(:)); dat(iX) = 0;
                [ixY(3) ixX(3)]=ind2sub(size(dat),iX);

                % find coor of 2by2 peak
                if round(mean(ixX)) == max(ixX), ixX(4) = min(ixX); elseif round(mean(ixX)) == min(ixX), ixX(4) = max(ixX); else, iserr2by2_='split peak'; end
                if round(mean(ixY)) == max(ixY), ixY(4) = min(ixY); elseif round(mean(ixY)) == min(ixY), ixY(4) = max(ixY); else, iserr2by2_='split peak'; end

                if ~isempty(iserr2by2_) % split peaks
                    dat=data3by3;
                    if ixY(2)== 2 && ixX(2) == 1 % cases 2 & 4
                        ixX(3:4) = [1,2];
                        if dat(1,1)+dat(1,2) > dat(3,1)+dat(3,2)
                            ixY(3:4) = 1; 
                        else
                            ixY(3:4) = 3; 
                        end
                    elseif ixY(2)== 2 && ixX(2) == 3 % cases 2 & 4
                        ixX(3:4) = [2,3];
                        if dat(1,2)+dat(1,3) > dat(3,2)+dat(3,3)
                            ixY(3:4) = 1; 
                        else
                            ixY(3:4) = 3; 
                        end                            
                    elseif ixX(2)== 2 && ixY(2) == 1 % cases 2 & 4
                        ixY(3:4) = [1,2];
                        if dat(1,1)+dat(2,1) > dat(1,3)+dat(2,3)
                            ixX(3:4) = 1; 
                        else
                            ixX(3:4) = 3; 
                        end
                    elseif ixX(2)== 2 && ixY(2) == 3 % cases 2 & 4
                        ixY(3:4) = [2,3];
                        if dat(2,1)+dat(3,1) > dat(2,3)+dat(3,3)
                            ixX(3:4) = 1; 
                        else
                            ixX(3:4) = 3; 
                        end
                    elseif abs(ixY(2)-2)+abs(ixX(2)-2) == 2  % cases 7 & 8
                        ixX(3:4) = fliplr(ixX(1:2));
                        ixY(3:4) = ixY(1:2);    
                    end
                    disp(iserr2by2_)
                    iserr2by2_ = [];
                end
                linIx = sort(sub2ind(size(dat),ixY,ixX)); % linear index

                [ytemp1 xtemp1] = ind2sub([3 3],linIx(1));
                ytemp = pxYc + ytemp1 - 2;
                xtemp = pxXc + xtemp1 - 2;
                temp = zeros(5); temp(ytemp,xtemp)=1; temp = padarray(temp,[2 2]);
                Window2by2 = Window5by5(ytemp:ytemp+1, xtemp:xtemp+1);
                [yCoor2by2(sp) xCoor2by2(sp)] = find(temp==1); % rectangle coord for 2 by 2
                % 2nd peak
                data3by3temp = data3by3;
                data3by3temp(linIx) = 0;                
                [v iX]=max(data3by3temp(:)); 
                [yPre2temp xPre2temp]=ind2sub(size(data3by3temp),iX);
                yPre2temp = pxYc + yPre2temp - 2;
                xPre2temp = pxXc + xPre2temp - 2;
                temp = zeros(5); temp(yPre2temp,xPre2temp)=1; temp = padarray(temp,[2 2]);
                [yPre2disp(sp) xPre2disp(sp)] = find(temp==1); % rectangle coord for 2 by 2
                
                % 2by2 localization
                is2by2localization = 0;
                if is2by2localization
                    for j = 1: size(XY,1)
                        fun2by2 = @(c,x) 0+c*exp(-intP1G*((x(:,1)-XY(j,1))/sg2by2/sqrt(2)).^2-intP1G*((x(:,2)-XY(j,2))/sg2by2/sqrt(2)).^2);
                        [I2by2(j),err2by2_(j),err2by2res_(:,j),EXITFLAG_] = gausFit(Window2by2,fun2by2,1); % ,c01G5by5,lb1G5by5,ub1G5by5);
                        %disp(EXITFLAG_)
                    end      
                    [v ix ]=min(err2by2_(:));
                    err2by2(sp) = err2by2_(ix);
                    %+ [pxXc pxYc] -1
                    PX2by2(sp) = XY(ix,1) + xCoor2by2(sp)-3;  
                    PY2by2(sp) = XY(ix,2) + yCoor2by2(sp)-3;  
                end
            end    
            %% gauss fit 2/2
            PX(sp)=nan;PY(sp)=nan; cc2_2(:,sp)=nan(11,1);
            if isTryFit(sp) % 2by2 peak
                tP2G_2 = inf;
                tP2G_1 = 1.75;
                tP2G = [tP2G_1 tP2G_1 tP2G_2 tP2G_2];
                xPre(sp) = xCoor2by2(sp)+0.5-2;
                yPre(sp) = yCoor2by2(sp)+0.5-2;
                xPre2(sp) = xPre2disp(sp)-2;
                yPre2(sp) = yPre2disp(sp)-2;
                xPre(sp) = 3; yPre(sp) = 3;
                %xPre2(sp) = 3; yPre2(sp) = 3;
                tol = inf;
                tol2 = inf; % background fit
                I2 = Window5by5(yPre2(sp),xPre2(sp)); 
                %if I2 < 0, tol2=1/tol; end;
                if I2 < 0, tol2=-inf; end;
                c02G5by5([3,5,8,10]) = [xPre(sp),yPre(sp),xPre2(sp),yPre2(sp)]; lb2G5by5([3,5,8,10]) = [xPre(sp),yPre(sp),xPre2(sp),yPre2(sp)]-tP2G; ub2G5by5([3,5,8,10]) = [xPre(sp),yPre(sp),xPre2(sp),yPre2(sp)]+tP2G;  
                c02G5by5([1,2,7]) = [0 1 0];  lb2G5by5([1,2,7]) = [0 0 I2/tol2];  ub2G5by5([1,2,7]) = [1*tol 1*tol I2*tol2];  % intensity
                %c02G5by5([3,5,8,10]) = [xPre(sp),yPre(sp),0,0]; lb2G5by5([3,5,8,10]) = [xPre(sp),yPre(sp),-inf,-inf]-tP2G; ub2G5by5([3,5,8,10]) = [xPre(sp),yPre(sp),inf,inf]+tP2G;  
                % FIT
                [cc2_2(:,sp),err2G5by5_2(sp),err2res_,EXITFLAG(3)] = gausFit(Window5by5,fun2_2,c02G5by5,lb2G5by5,ub2G5by5);
                
                % ERR
                %err2G5by5_2(sp) = (sqrt(err2G5by5_2(sp))/numel(Window5by5) ).^(1/intP2G)*100;
                
                %% the center
                PX_ = cc2_2(3,sp);
                PY_ = cc2_2(5,sp);
                PX(sp) = PX_;
                PY(sp) = PY_;     
            end

            %% calc XY value on the image (xyzDataGaus)
            Xgaus(iyTrack(i),ixTrack(i)) =  double(xBW(iyTrack(i),ixTrack(i))) + PX(sp) -px0;
            Ygaus(iyTrack(i),ixTrack(i)) =  double(yBW(iyTrack(i),ixTrack(i))) + PY(sp) -px0;
            Zgaus(iyTrack(i),ixTrack(i)) =  zBW(iyTrack(i),ixTrack(i));
            %INTmax(iyTrack(i),ixTrack(i)) = intMax(sp); % max window
            %INTpeak(iyTrack(i),ixTrack(i)) = intPeak_(sp); % max window - bckgrnd
            %INTgaus(iyTrack(i),ixTrack(i)) = cc2_2(2,sp); % gausfit to
            %normalized window
            INTnorm(iyTrack(i),ixTrack(i)) = intMax(sp)*intPeak_(sp)*cc2_2(2,sp); % actual intensity
            disp(sprintf('fr#:%i, sp#:%i, window center:%.02f,%.02f, EXITFLAG: %i ',ixTrack(i),sp,xBW(iyTrack(i),ixTrack(i)),yBW(iyTrack(i),ixTrack(i)),EXITFLAG(3)));

            %% for display only
            viewSpot2dispLoop;

            waitbar(i/nSpots);
            sp2(sp) = i;
            sp = sp + 1;
        end
        SP = sp-1;
        %save(xyzDataGausFileNm,'Xgaus','Ygaus','Xgaus2by2','Ygaus2by2','X_','Y_');
        save(xyzDataGausFileNm,'Xgaus','Ygaus','Zgaus','INTnorm');
        viewSpot3dispSave;
        close(h);
    end
    
    if isGausFitLocalize && ~exist(xyzDataGausFiltFileNm)
        assignFileNames
        load(xyzDataGausFileNm)
        load(spotWinFileNm, 'xBW','yBW')
        % define output XY array
        [iyTrack ixTrack] = find(xBW>0); % positions of spot data in the X array
        
        % clear out nonfit spots 
        if ~exist('spotSelVal.mat') % by selecting the threshold values
            viewSpot4dispPlot; % outputs=>spotSelVal.mat & ixSel (loads from viewSpot3dispSave and plots)
            return;
            load(xyzDataGausFileNm)
        else        % using the saved threshold values
            load spotSelVal % 'ixSel','ixNonSel','errThreshold','intThreshold'
            load spotInfo
            ixInt5 = find(intFit>=intThreshold); % selected
            ixInt5_2 = find(intFit<intThreshold);
            err5 = err2G5by5_2(ixInt5_2);
            ixErr5 = find(err5<=errThreshold); % also selected
            ixSel = sort([ixInt5_2(ixErr5) ixInt5]); % passing spots
        end
        disp(sprintf('GAUSS FIT elimination: %.02f percent of the spots are discarded',(1-(numel(ixSel)/numel(intFit)))*100));
        fprintf(fidLog,'GAUSS FIT elimination: %.02f percent of the spots are discarded',(1-(numel(ixSel)/numel(intFit)))*100)
        
        % filter
        gausFitSel  = zeros(size(Xgaus));
        gausFitSel(iyTrack(sp2(ixSel)),ixTrack(sp2(ixSel))) = 1;
        Xgaus = Xgaus.*gausFitSel;
        Ygaus = Ygaus.*gausFitSel;
        Zgaus = Zgaus.*gausFitSel;
                
        X = Xgaus; Y = Ygaus; Z = Zgaus; clear Xgaus Ygaus;
        INT = INTnorm; clear INTnorm;
        save(xyzDataGausFiltFileNm,'X','Y','Z','INT');
        delete(traceDataFileNm0)
    end
    
    
    %% TIME To FIND OUT THE TRACES
    if ~exist(traceDataFileNm0)
        load(posDataFileNm);
        if exist(xyzDataGausFiltFileNm)
            load(xyzDataGausFiltFileNm);
        elseif exist(xyzDataGausFileNm)
            error('GausFiltFile missing')
            load(xyzDataGausFileNm); 
        else
            error('GausFiltFile missing')
            load(xyzDataFileNm); 
        end
        
        
        if ~exist('Z'), Z = zeros(size(X)); end;
        [iy ix] = find(X>0); debugSpt = [X(X>0) Y(Y>0) ix];
        
        isFindRecruitment=1;
        if isFindRecruitment
            sptJmp = 1;
            sptReAppearTime = 2; 
            minTraceLength = 2;
            %INT = zeros(size(X));
        end        
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
                        for n=1:Boy % all spots
                            dif(n,l)=sqrt((AslX-X(n,k+l-1))^2+(AslY-Y(n,k+l-1))^2+[(AslZ-Z(n,k+l-1))*PlaneDist/PixelSize]^2);        
                        end
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
        save(traceDataFileNm0,'TraceX','TraceY','TraceZ','TraceINT'); 
        save(cfgTraceFileNm,'sptJmp','sptReAppearTime','minTraceLength');
        delete(traceDataFileNm)
    end
    
    
    if 0 && ~exist(traceDataFileNm) && isCombTraces
        %% COMBINE TRACES
        load(traceDataFileNm0)
        [m n]=size(TraceX);
        h = waitbar(0,'Combining the traces...');
        disp('Combining the traces...');

        [Boy2,Frames]=size(TraceX);
        tic;
        com=0;
        for i=1:Boy2-1
        %     TraceDif1=0;TraceDif2=0;TraceDif3=0;
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
                        if min(TraceDif)<10
                            for k=1:Frames
                                if [TraceINT(i,k)+TraceINT(j,k)] > 0
                                TraceX(i,k)= [TraceX(i,k)*TraceINT(i,k)+TraceX(j,k)*TraceINT(j,k)]/[TraceINT(i,k)+TraceINT(j,k)];
                                TraceY(i,k)= [TraceY(i,k)*TraceINT(i,k)+TraceY(j,k)*TraceINT(j,k)]/[TraceINT(i,k)+TraceINT(j,k)];
                                TraceZ(i,k)= [TraceZ(i,k)*TraceINT(i,k)+TraceZ(j,k)*TraceINT(j,k)]/[TraceINT(i,k)+TraceINT(j,k)];
                                TraceINT(i,k)= [TraceINT(i,k)*TraceINT(i,k)+TraceINT(j,k)*TraceINT(j,k)]/[TraceINT(i,k)+TraceINT(j,k)];
                    %             TnaceX(i,k)=TraceX(i,k); TnaceY(i,k)=TraceY(i,k); TnaceX(j,k)=TraceX(j,k); TraceY(j,k)=TraceY(j,k);
                                end
                            end
                            TraceX(j,:)=NaN; TraceY(j,:)=NaN; TraceZ(j,:)=NaN; TraceINT(j,:)=NaN;
                            LastElement=max(find(TraceX(i,:)>0));
                            LastBefore=LastElement-1;
                            com=com+1;
                        end
                    end
                end
            end
                %Tcomp0 = toc;
        end
        waitbar(i / Boy2)
        close(h)
        % 'traceData-coeff%d.mat'
        Tcomp = toc;
        save(traceDataFileNm,'TraceX','TraceY','TraceZ','TraceINT','Tcomp')
        
    end
    

% PLOT TRACES
    disp('preparing plot data...')
    if isCombTraces
        load(traceDataFileNm);
    else
        load(traceDataFileNm0);
    end
        if exist(xyzDataGausFiltFileNm)
            load(xyzDataGausFiltFileNm);
        elseif exist(xyzDataGausFileNm)
            load(xyzDataGausFileNm); 
        else
            load(xyzDataFileNm); 
        end
    
    for tr = 1:size(TraceX,1)
        %% fill gaps (missing frames in jumpy traces)
        traceX = TraceX(tr,:);
        traceY = TraceY(tr,:);
        traceZ = TraceZ(tr,:);
        ind = find(traceX~=0);
        frst = min(ind); last = max(ind);
        %iJump = find(trace(frst:last)==0)+frst-1;

        jmp = abs((traceX(frst:last)>0)-1);
        if sum(jmp > 0)
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
    end
    
    isCropData = 0; 
    if isCropData
        Tr1=1;
        Tr2=300;
        %Tr2=size(TraceY,1);
        tt1 = 1; 
        tt2 = size(TraceY,2);
        tt2 = 50;
        dTr = 1;
        TraceX = TraceX(Tr1:dTr:Tr2,tt1:tt2);
        TraceY = TraceY(Tr1:dTr:Tr2,tt1:tt2);
        TraceZ = TraceZ(Tr1:dTr:Tr2,tt1:tt2);
    end
        
    % color data
    CplotVecN = size(TraceX,1); % # of traces
    Nframe = size(TraceX,2); % # of frames
    useCData = 1;
    CData = TraceZ;
    CData = 1:size(TraceX,2); CData=repmat(CData,[size(TraceX,1) 1]);
    if useCData 
        Zmin = min(CData(CData~=0));
        Zmax = max(CData(:));
        Zrange = Zmax - Zmin;
        Cplot = round(63*(CData-Zmin) / Zrange)+1;
        CM = colormap('jet');
    else
        CplotVec = 1: CplotVecN;
        CplotVec = rem(CplotVec,64)+1;
        Cplot = repmat(CplotVec',1,Nframe);
        CM = colormap('lines');
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
    padSize = 0;
    if exist(xyzDataGausFileNm), load(xyzDataGausFileNm); else, load(xyzDataFileNm); end;
    zeroIx = find(X == 0);
    X = X + padSize; Y = Y + padSize;
    X(zeroIx) = 0; Y(zeroIx) = 0;
    zeroIx = find(TraceX == 0);
    TraceX = TraceX + padSize; TraceY = TraceY + padSize;
    TraceX(zeroIx) = 0; TraceY(zeroIx) = 0;
    En = Frames;
    [Boy2]=size(TraceX,1);
    
    nonZeroIx = find(TraceY>0);
    dispPx = {'*',2}; % plot3k
    minimum = 1;
    maximum = max(Cplot(:));

%% PLOT    
    NAN = find(isnan(TraceX));
    TraceX(NAN)=0;
    TraceY(NAN)=0;
    TraceZ(NAN)=0;
    

    % parameters:
    load(posDataFileNm,'sptReAppearTime'); %(frames) use the value from tracker function generating trace values

    load fname;
    imgZFout = 'TraceImage';

    img2D = imread(fname,1); 
    img2D = img2D(yy1:yy2,xx1:xx2);
    img2D = padarray(img2D,[padSize padSize]);
    imgFrm = uint16(zeros(size(img2D))); % image padding frame
    imgFrm(yy1+3:yy2-3,xx1+3:xx2-3)=1;    
    imSz = size(img2D');
    TraceY(nonZeroIx) = imSz(2)-TraceY(nonZeroIx)+1;

    mag = 5;
    [mag, pos, m, n ] = calcMaxMag(img2D,mag);
    colormap('gray');

    zMax = max(TraceZ(:));
    if zMax == 0, zMax = 1; TraceZ = TraceZ+1; end;

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
    lastZ = TraceZ(:,1);
    frm1 = 40;
    frm1 = 1;
    tit = 'image';
    m_ = uint16(imSz(1)); n_ = uint16(imSz(2));
    %pos=get(0,'ScreenSize');
    %pos=uint16(pos(3:4)) - [m n-35];
    fig=figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos/2 m n]);
    axe=axes('Parent',fig,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren','XLim',0.5+[0 m],'YLim',0.5+[0 n]);
    hWB =  waitbar(0,'marking spots...');
    for ixFrm = frm1:En
        img2D = imread(fname,ixFrm); 
        img2D = flipud(img2D(yy1:yy2,xx1:xx2));
        % padding
        img2D = padarray(img2D,[padSize padSize]);
        img2D = img2D.*imgFrm;
        %axe=axes('Parent',fig,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','XLim',0.5+[0 m],'YLim',0.5+[0 n]);
        hold on
        set(axe,'Parent',fig,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','XLim',1+[0 m_],'YLim',1+[0 n_]);
        hImg = imagesc(img2D,'Parent',axe); %axis image; 
        axis tight
        ixTrc = find(TraceX(:,ixFrm)>0);
        currX = TraceX(:,ixFrm);
        currY = TraceY(:,ixFrm);
        currZ = TraceZ(:,ixFrm);
        uistack(hImg,'bottom');
        hold on
        ixTrc2 = find(lastX(ixTrc)>0);
        ix = ixTrc(ixTrc2); % indices of traces tracked in the current frame    
        dspTrcIx = find(ixFrm == dspTrcFrm); % index for dissappearing traces
        showTrace =1; % puts arrows
        isShowTrace = 1;
        if showTrace && isShowTrace
            hQdel = hQ(dspTrcIx,:); % handles for discont. traces
            delete(hQdel(find(hQdel~=0)));    % remove the traces of the discontinued traces
        end
        isColorbyTime = 1;
        if isShowTrace
            for i = 1:round(length(ix)) % for each trace 
                iL = ix(i);  % index for each line
                hQ(iL,ixFrm-frm1+1)=quiver(lastX(iL),lastY(iL),currX(iL)-lastX(iL),currY(iL)-lastY(iL),'Color',CM(Cplot(iL,ixFrm),:));
                if ~showTrace
                    adjust_quiver_arrowhead_size(hQ(iL,ixFrm-frm1+1),5)
                end
            end
        end
        hscat = scatter(X(find(Y(:,ixFrm)>0),ixFrm),imSz(2)-Y(find(Y(:,ixFrm)>0),ixFrm)+1,50,'r','d','LineWidth',0.5);

        % print images
        imgFig = getframe(gcf);
        imgOut = imgFig.cdata;
        figPos = get(gcf,'Position');
        %imgOut = imgOut((figPos(4)-m)/2+1:(figPos(4)+m)/2,(figPos(3)-n)/2+1:(figPos(3)+n)/2,:);
        imgOut = imgOut((figPos(4)-n)/2+1:(figPos(4)+n)/2,(figPos(3)-m)/2+1:(figPos(3)+m)/2,:);
        
        fOut = sprintf('%s.tif',imgZFout);
        if ixFrm == 2
            if exist([imgZFout '.tif'])
                delete([imgZFout '.tif']);
            end
            imwrite(imgOut,fOut,'Compression', 'none') 
        else
            imwrite(imgOut,fOut,'WriteMode','append','Compression', 'none') 
        end
        delete(hImg)
        delete(hscat)

        %waitforbuttonpress
        if ~showTrace
            delete(hQ(find(hQ(:,ixFrm-1)~=0),ixFrm-1));
        end

        %export_fig(imgZFout,'-append');
        %fprintf('frame %i/%i \n',ixFrm,En);
        lastX = currX; 
        lastY = currY;
        lastZ = currZ;
        waitbar(ixFrm/En)
    end
    close(hWB);
    hold off; % release the figure 
    
    return;
    
    %tiff2stack('TraceImage'); % combine frames in a stack
    isTrace = 1;
    save('isTraceFile','isTrace');
    recruitmentTrack;
    isTrace = 0;
    save('isTraceFile','isTrace');    
    recruitmentTrack;
    delete('isTraceFile.mat');
    
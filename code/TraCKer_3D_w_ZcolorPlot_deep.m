
% 'style','text','BackgroundColor',[1 1 1],
    clear; close all;
    F = findall(0,'type','figure','tag','TMWWaitbar'); delete(F)
    isCropXY = 0;
    isDebugFilt = 0; 
    isRangeLow = 0;
    noTile = 1;
    binFrameTime = 0; % [ms]
    isCombTraces = 0;
    maxTraceSpeed = nan; % [um/sec] used for calculation of allowed trace jump(not used)
%    minXYspread = 15; % [px] std of the xy location along the trace
    minXYspread = 0.8; % [px] std of the xy location along the trace
    minXYspread = nan;
    
%minXYspread = 1.6; % [px] std of the xy location along the trace
    
    
    
    % compressed sensing parameters
    div = 2; % CS "zoom-factor"
    eps = 1.5; % optimization multiplicative constaint parameter
    photon_per_count = 1/300; % := 1/gain
    ccd_base_line = 2000; 
    
    
    inFocus = 1; % for detection gaussian size. if sharp spots:1, ow:0
    inFocus_w_1X = 0; % if 1.5x mag not used 
    isDiffData = 0;

    cfg = struct;
    cfg.report = struct;
    
    isBALM = 0;
    isOptimFit = 0;
    isGausFit = 0; % Gaussian localization
    %isGausFitLocalize = 1; % gaussian fit filtering and localization
    StackNum = 1;   
    
    if inFocus_w_1X | isDiffData | inFocus
        WindowSize = 3; 
        BigWindowSize=WindowSize+2;
    else
        WindowSize = 5; 
        %BigWindowSize=WindowSize+4;
        BigWindowSize=WindowSize;
    end
    
    coeffFitParamCoverage = 1000; % number of frames for a fit param
%coeffFitParamCoverage = 50; % number of frames for a fit param
    
        
    %% select frames
    frstFrm = 1; % cropping
    frstFrm2 = 1; 
    lastFrm = [];
    if exist('frames.txt') % frames already selected
        fo = fopen('frames.txt');
        inp = fgetl(fo);
        sc = strfind(inp,':'); % pos. semicolon
        frm1 = str2num(inp(1:sc-1));
        frm2 = str2num(inp(sc+1:end));
        frstFrm2 = frm1; % sugg. trimming the beginning of the movie
        lastFrm = frm2;
        disp(sprintf('selected frames: frames %i-%i',frm1,frm2));
        fclose(fo);
    elseif exist('BWselect.mat') % background selected
        
        isFrmCrp = input('display intensity change?[y/n]','s');
        if strcmp(isFrmCrp,'y')
            selFrames
        else
            return;
        end
    end
    
    %% read info files
    % read cell info
    infoFolder = brdir('info');
    cfg.cell = struct;
    fcell = [infoFolder '/info-cell.txt'];
    if isempty(fcell), error('no cell info is found'); end;
    hfcell = fopen(fcell);
    cfg.cell.cellType = fscanf(hfcell,'cellType:%s\n');
    cfg.cell.freezeDate = fscanf(hfcell,'freezeDate:%s\n'); % eg 10h00m150210
    cfg.cell.thawDate = fscanf(hfcell,'thawDate:%s\n'); % eg 10h00m150210
    cfg.cell.platingTime = fscanf(hfcell,'platingTime:%s\n'); % eg 10h00m150210
    cfg.cell.incubOutTime = fscanf(hfcell,'incubOutTime:%s\n'); % eg 10h00m150210
    cfg.cell.surfTreat = fscanf(hfcell,'surfTreat:%s\n');
    cfg.cell.chemTreat = fscanf(hfcell,'chemTreat:%s\n');
    cfg.cell.incubMed = fscanf(hfcell,'incubMed:%s\n');
    cfg.cell.imageMed = fscanf(hfcell,'microsMed:%s\n');
    cfg.cell.confluency = fscanf(hfcell,'confluency:%s\n');
    cfg.cell.maxConfluency = fscanf(hfcell,'maxConfluency:%s\n');
    cfg.cell.other = fscanf(hfcell,'other:%s');
    fclose(hfcell);
    
    % read exp info
    cfg.exp = struct;
    fexp = [infoFolder '/info-exp.txt'];
    if isempty(fexp), error('no exp info is found'); end;
    hfexp = fopen(fexp);
    cfg.exp.date = fscanf(hfexp,'date:%s\n');
    cfg.exp.time = fscanf(hfexp,'time:%s\n');
    cfg.exp.expSetup = fscanf(hfexp,'expSetup:%s\n');
    cfg.exp.other = fscanf(hfexp,'other:%s');
    fclose(hfexp);
    
        
    % read acquisition info
    cfg.acq = struct;
    facq = [infoFolder '/info-acq.txt'];
    if isempty(facq), error('no acquisition info is found'); end;
    hfacq = fopen(facq);
    cfg.acq.ROIcam = fscanf(hfacq,'camROI: %s\n');
    cfg.acq.ROInd = fscanf(hfacq,'ndROI:%s\n');
    cfg.acq.TIRFstagePosEPI = fscanf(hfacq,'EPIstagePos:%s\n');
    cfg.acq.TIRFstagePosTIRF = fscanf(hfacq,'TIRFstagePos:%s\n');
    cfg.acq.tempBeforePositioning = fscanf(hfacq,'1-tempBeforePositioning:%s\n');
    cfg.acq.tempChangeInTenMin = fscanf(hfacq,'2-tempChangeInTenMin:%s\n');
    cfg.acq.preImg = fscanf(hfacq,'preImg:%s\n'); % mJoule
    cfg.acq.adjDur = fscanf(hfacq,'calib:%s\n'); % [sec] z and tirf angle adjustment (needs to be taken out)
    cfg.acq.acqTime = fscanf(hfacq,'acqTime:%s\n'); % [ms]
    cfg.acq.power = fscanf(hfacq,'power:%s\n');
    cfg.acq.gain = fscanf(hfacq,'gain:%s\n');
    cfg.acq.NDfilt = fscanf(hfacq,'NDfilt:%s\n');
    cfg.acq.pixelCAM = fscanf(hfacq,'pixelCAM:%s\n');
    cfg.acq.magObj = fscanf(hfacq,'magObj:%s\n');
    cfg.acq.magEx	 = fscanf(hfacq,'ExMag:%s\n');
    cfg.acq.other = fscanf(hfacq,'other:%s');
    fclose(hfacq);
    
    cfg.proc.frm1=frm1;
    
    
    if 0
        % read process info
        cfg.proc = struct;
        fproc = [infoFolder '/info-proc.txt'];
        if isempty(fproc), error('no bleach info is found'); end;
        hfproc = fopen(fproc);
        cfg.proc.trim = fscanf(hfproc,'trim:%s\n'); 
        fclose(hfproc);

        frstFrm2 =    str2num(cfg.proc.trim);
        if ~frstFrm2
            frstFrm2 = 1;
        end
    end
    
    %% READ the ORIGINAL file
    if exist('fname.mat')
        load fname
        if ~exist('fname')
            fname = sprintf('../\%s',fname);
        end
    end
    
    %defineMask;
    
    acqTime = str2num(cfg.acq.acqTime)/1000; % sec
    if exist('acqTime.txt'), acqFn=fopen('acqTime.txt'); acqTime = str2num(fscanf(acqFn,'%s'))/1000; end;


    % upd binFrame (binning of the frames)
    if binFrameTime
        binFrame = round(binFrameTime/acqTime/1000); % [frames]
    else
        binFrame = 1; % [frames] no binnning
    end
    fnameBck = [fname(1:end-4) '_background.tif']; % default background file
    label = fname(end-6:end);
    if ~exist(fnameBck)
        fTif1 = dir('../*.tif');
        fTif2 = dir('../../*.tif');
        fTif = [fTif1 ;fTif2];
        for i = 1:length(fTif)
            fTif(i).name;
            if ~isempty(strfind(fTif(i).name,label)) && isempty(strfind(fTif(i).name,'MAX'))
                fnameBck = ['../' fTif(i).name];
                if ~exist(fnameBck)
                    fnameBck = ['../../' fTif(i).name];
                end
            end
        end
        if ~exist(fnameBck)
            fnameBck = fname;
        end
    else % no cropping 
        imgTemp = imread(fnameBck,1);
        BWselect = ones(size(imgTemp));
        save('BWselect','BWselect');
    end
    
    posData_File =dir('posData-coeff*');
    
    %% data filename and image info
    if ~exist('fname')
        ch1_1=strfind(fname,'\'); ch1_2=strfind(fname,'/');
        if length(ch1_1) > length(ch1_2), ch1 = ch1_1; else ch1 = ch1_2; end;
        ch2=strfind(fname,'x');
        ch3=strfind(fname,'xy');
        En1= str2num(fname(ch1(end)+1:ch2(1)-1));
        Boy1= str2num(fname(ch2(1)+1:ch3(1)-1));
    elseif ~isempty(posData_File)
        imgFrst = imread(fname);
        [Boy1,En1]=size(imgFrst);  
        load(posData_File.name,'Frames');
    else
        imgFrst = imread(fname);
        [Boy1,En1]=size(imgFrst);    
        dirFN = dir(fname);
        if dirFN.bytes > 1e10
            Frames=findNumFrames(fname);
        else
            infFN = imfinfo(fname);
            Frames = numel(infFN);
        end
        clear infFN;
        
        Frames = floor(Frames/binFrame);
        %Frames = 5500;
        % Frames = 100; % # of frames to be analyzed
        if ~isempty(lastFrm)
            Frames = lastFrm;
        else
            lastFrm = Frames;
        end
    end
    
    if isCropXY
        xx1 = 30; xx2 = 65; % crop
        yy1 = 20; yy2 = 65; % crop
        En1 = xx2-xx1+1;
        Boy1 = yy2-yy1+1;
    else
        xx1 = 1; xx2 = En1; % crop
        yy1 = 1; yy2 = Boy1; % crop
    end    
    sizeImg=[Boy1,En1];
    
    % crop number
    CD=cd;
    if isunix,  slashPos= strfind(CD,'/'); else slashPos= strfind(CD,'\'); end;
    folderName=CD(slashPos(end)+1:end);
    cropNo=sscanf(folderName,'crop%i');
    
    cfg.img = struct;
    cfg.img.boy = Boy1;
    cfg.img.en = En1;
    cfg.img.fname = fname;
    cfg.img.frstFrm = frstFrm2;
    cfg.img.lastFrm = lastFrm;
    cfg.img.binFrame = binFrame;
    cfg.img.binFrameTime_ms = binFrameTime;    
    cfg.fit = struct;

    
    cfg.cs = struct; % compressed sensing
    cfg.cs.div = div;
    cfg.cs.eps = eps;
    cfg.cs.photon_per_count = photon_per_count;
    cfg.cs.ccd_base_line = ccd_base_line;
    
    coeffMat = 'CoeffFit.mat';
    Coeff = 0;
    CoeffFit = 0;
    
    coeffFound = 0;
    if exist(coeffMat)
        load(coeffMat);
        
        coeffFound = 1;
    else
        display('No Coeff file was found');
    end
    assignFileNames

    %% FIND X & Y
    if ~exist(xyzDataGausFileNm) 
       
        %% detection parameters
        topHat4by4 = [-0.0833   -0.0833   -0.0833   -0.0833;-0.0833    0.2500    0.2500   -0.0833;-0.0833    0.2500    0.2500   -0.0833;-0.0833   -0.0833   -0.0833   -0.0833]; % topHat
        topHat5by5 = [-0.0625   -0.0625   -0.0625   -0.0625   -0.0625;-0.0625    0.1111    0.1111    0.1111   -0.0625;-0.0625    0.1111    0.1111    0.1111   -0.0625;-0.0625    0.1111    0.1111    0.1111   -0.0625;-0.0625   -0.0625   -0.0625   -0.0625   -0.0625];
        topHat5by5 = [-0.0649   -0.0649   -0.0649   -0.0649   -0.0649; -0.0649    0.1087    0.1087    0.1087   -0.0649;-0.0649    0.1087    0.1087    0.1087   -0.0649;-0.0649    0.1087    0.1087    0.1087   -0.0649;-0.0649   -0.0649   -0.0649   -0.0649   -0.0649];
        topHat7by7 = [-0.0417   -0.0417   -0.0417   -0.0417   -0.0417   -0.0417   -0.0417;-0.0417    0.0400    0.0400    0.0400    0.0400    0.0400   -0.0417;-0.0417    0.0400    0.0400    0.0400    0.0400    0.0400   -0.0417;-0.0417    0.0400    0.0400    0.0400    0.0400    0.0400   -0.0417;-0.0417    0.0400    0.0400    0.0400    0.0400    0.0400   -0.0417;-0.0417    0.0400    0.0400    0.0400    0.0400    0.0400   -0.0417;-0.0417   -0.0417   -0.0417   -0.0417   -0.0417   -0.0417   -0.0417];
        lap=[-1,-1,-1;-1,8,-1;-1,-1,-1];
        lap5x5= -1/16*ones(5);
        lap5x5(2:4,2:4)=0;
        lap5x5(3,3)=1;
        lap5x5 = lap5x5*8;
        sqr2x2= ones(2)/4;

        isOld = 0;
        
        if inFocus_w_1X
            gausKernelSz = 3;
            gausKernelSg = .5;
            gausKernelSzPadded = 3;
            %lap = [-1    -1    -1    -1;-1     3     3    -1; -1     3     3    -1; -1    -1    -1    -1];
            %lap =1;
        elseif inFocus
            gausKernelSz = 3;
            gausKernelSg = 0.7;
            gausKernelSzPadded = 5;
        elseif isOld
            gausKernelSz = 3;
            gausKernelSg = 0.7;
            gausKernelSzPadded = 7;
        else
            gausKernelSz = 5;
            gausKernelSg = 1;
            gausKernelSzPadded = 7;
        end
%         gausKernelSzPadded = 9;
%         gausKernelSz = 5;
%         gausKernelSg = 1;
        
        pdSz = (gausKernelSzPadded - gausKernelSz)/2; % pad size
        gausHat1=fspecial('gaussian', gausKernelSz, gausKernelSg);
        cfg.fit.detect = struct;
        cfg.acq.acqTime = acqTime; % sec
        
        elev=1/(gausKernelSzPadded^2-gausKernelSz^2);        
        if elev==inf, elev=0; end;
        gausHat1 = gausHat1 + elev;
        gausHatPad1 = padarray(gausHat1,[pdSz pdSz]);
        gausHatPad1 = gausHatPad1 - elev;
        gausHat1 = gausHatPad1;
        gausHat2 = gausHat1;
        
        cfg.fit.detect.gausKernelSzPadded = gausKernelSzPadded;
        cfg.fit.detect.WindowSize = WindowSize;
        cfg.fit.detect.BigWindowSize = BigWindowSize;
        cfg.fit.detect.gausKernelSg = gausKernelSg;
        cfg.fit.detect.gausKernelSz = gausKernelSz;
        cfg.fit.detect.gausHat1 = gausHat1;
        
        %% set the coefficient
        skipAdjCoeff = 0;
        if coeffFound
            skipAdjCoeff = 1;
            if 0
                ans=input('do you want to adjust the coeff? (y/n)','s')
                if strcmp(ans,'y')
                    skipAdjCoeff = 0;
                elseif strcmp(ans,'n')
                    skipAdjCoeff = 1;
                else
                    return;
                end
            end
        end
            
        %if ~coeffFound && ~isOptimFit
        if ~skipAdjCoeff && ~isOptimFit
            %% use coeff in the later crops
            CoeffSugg=[];
            if cropNo > 1 % if multi crop files of the data set the Coeff for the following crops
                cd(sprintf('%scrop%i',CD(1:slashPos(end)),cropNo-1))
                load CoeffFit
                CoeffSugg = CoeffFit(end); % suggested coeffiecient
                clear CoeffFit; 
                cd(CD)
            end
            if coeffFound % set parameters acc. to the loaded data
                load CoeffFit
                CoeffAuto = CoeffFit/max(CoeffFit);
                CoeffLoaded = CoeffFit(1);
            else
                CoeffAuto = ones(Frames-frstFrm2+1,1);
            end
                       
            %% generate preview image
            nTile = 1;
            CoeffFit = CoeffAuto;
            genPreviewImg;  % generate the image SHOW
            
            %% image figure
            imSize = size(SHOW);
            m = imSize(1); n = imSize(2);
            tit = 'set the coefficient values for peak detection';
            hFigDisp=figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256));
            set(gcf,'units','normalized','outerposition',[0 0 1 1])
            kk =0.1;
            axe=axes('Parent',hFigDisp,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','outerposition',[0 0 1 1],'Nextplot','replacechildren','XLim',0.5+[0 m],'YLim',0.5+[0 n]);

            %mx = max([PREdataFiltTiled1(:); PREdataFiltTiled2(:)]);
            %Scale=[1:1:mx];
            
            %% control figure
            tit2 = 'controls';
            hFigCont = figure('DoubleBuffer','on','Menubar','none','Name',tit2,'NumberTitle','off','Colormap',gray(256));
            figPos = [1921 -461 1080 1844];
            %figPos = [1920 -245 1080 1844];
            set(hFigCont,'Position', figPos);

            % subplot(1/2) draw coeff norm factor
            hCoeffPlot = subplot(2,1,1);
            coeffPlot = plot(CoeffAuto); title('Coefficient change'); ylabel('Coefficient Norm.'); xlabel('time');
            
            % coeff parameter settings
            CoeffFitMax=1.2;
            nFrames = Frames-frstFrm2+1;
            numCoeffFitParam = floor(nFrames/coeffFitParamCoverage); % number of coefficient fit parameters
            if ~numCoeffFitParam, numCoeffFitParam=1;end;
            axis tight;
            %ylim([0 CoeffFitMax]);
            set(hCoeffPlot,'Units','pixels');
            figPos = get(hCoeffPlot,'Position');
            xs = coeffFitParamCoverage/nFrames*(figPos(3));
            qCoeffSel = 0;
            yp = figPos(2);
            wp = 5;
            hp = figPos(4);
            hImgCoeffParamText0_ = uicontrol('style','text','BackgroundColor',[1 1 1],'String','Coeff:','Position',[figPos(1)-100 yp-40 60 15]);
            hImgCoeffParamText0 = uicontrol('style','text','BackgroundColor',[1 1 1],'String','0','Position',[figPos(1)-30 yp-40 60 15]);
            hBtnFrameSet0 = uicontrol('style','pushbutton','units','pixel','position',[figPos(1)-50 yp-65 80 20],'Callback', sprintf('activeSliderNo=0; toggleOut=1; nTile=%i; isNewTileImage = 1;',1),'String',sprintf('fr#:1')); 
            for i = 1:numCoeffFitParam-1
                frameSet = coeffFitParamCoverage*i;
                xp = figPos(1)+i*xs;
                hImgCoeffParam(i) = uicontrol('style','slider','units','pixel','position',[xp yp wp hp],'Value',CoeffAuto(frameSet),'Callback', sprintf('activeSliderNo=%i; isNewTileImage=1;toggleOut=1; nTile=%i;',i,frameSet));
                set(hImgCoeffParam(i),'Max',CoeffFitMax)
                hImgCoeffParamText(i) = uicontrol('style','text','BackgroundColor',[1 1 1],'String','0','Position',[xp-30 yp-40 60 15]);
                hBtnFrameSet(i) = uicontrol('style','pushbutton','units','pixel','position',[xp-50 yp-65 80 20],'Callback', sprintf('activeSliderNo=%i; toggleOut=1; nTile=%i; isNewTileImage = 1;',i,frameSet),'String',sprintf('fr#:%i',frameSet)); 
            end
            if isempty(i), i = 0;xp = figPos(1); end
            frameSet = Frames-frstFrm2+1;
            hImgCoeffParam(i+1) = uicontrol('style','slider','units','pixel','position',[figPos(1)+figPos(3) yp wp hp],'Value',CoeffAuto(frameSet),'Callback', sprintf('activeSliderNo=%i; isNewTileImage=1;toggleOut=1; nTile=%i;',i+1,frameSet)); 
            set(hImgCoeffParam(i+1),'Max',CoeffFitMax)
            hImgCoeffParamText(i+1) = uicontrol('style','text','BackgroundColor',[1 1 1],'String','0','Position',[figPos(1)+figPos(3)-30 yp-40 60 15]);
            hBtnFrameSet(i+1) = uicontrol('style','pushbutton','units','pixel','position',[figPos(1)+figPos(3)-50 yp-65 80 20],'Callback', sprintf('activeSliderNo=%i;toggleOut=1; nTile=%i; isNewTileImage = 1;',i+1,frameSet),'String',sprintf('fr#:%i',frameSet)); 

            % subplot(2/2) draw histogram
            axe=axes('Parent',hFigCont ,'Visible','off');
            for i = 1:50
                histImg(:,:,i) = imread(fname,i);
            end
            histImgMax = max(histImg,3);
            subplot(2,1,2)
            hist([0 ;double(histImgMax(histImgMax>0))],1000); histSHOW = hist([0; double(histImgMax(histImgMax>0))],1000);
            [ ~, histPeakPos] = max(histSHOW); ytick = get(gca,'YTick'); maxHistDisp = ytick(end);

            %% lookup table and initial coeff settings
            bt = yp-220;
            hTextCoeff = uicontrol('style','text','BackgroundColor',[1 1 1],'String','coefficient','Position',[20 bt+20 160 15]);
            hCoeff = uicontrol('style','slider','units','pixel','position',[20 bt 300 20],'Callback','activeSliderNo=0;'); 
            hTextImgMin = uicontrol('style','text','BackgroundColor',[1 1 1],'String','min','Position',[20 bt+60 160 15]);
            hImgMin = uicontrol('style','slider','units','pixel','position',[20 bt+40 300 20]); 
            hTextImgMax = uicontrol('style','text','BackgroundColor',[1 1 1],'String','max','Position',[20 bt+100 160 15]);
            hImgMax = uicontrol('style','slider','units','pixel','position',[20 bt+80 300 20]); 
            hTextImgGamma = uicontrol('style','text','BackgroundColor',[1 1 1],'String','gamma','Position',[20 bt+140 160 15]);
            hImgGamma = uicontrol('style','slider','units','pixel','position',[20 bt+120 300 20]); 
            hBtnCoeff = uicontrol('style','pushbutton','units','pixel','position',[340 bt 100 20],'Callback', 'q=1;','String','QUIT'); 
            hToggleBtnImgWin = uicontrol('style','togglebutton','units','pixel','position',[400 bt+40 100 100],'Callback', 'toggleOut=1;','String','browse mode'); 
            
            %hTextBckCoeff = uicontrol('style','text','String','Select coeff parameters','Position',[xp+400 yp-60 100 20]);
            % = uicontrol('style','pushbutton','units','pixel','position',[xp+400 yp-100 100 20],'Callback', 'qCoeffParamSel=1;','String','OK'); 
            
            %% initialize lookup table and coeff values
            imgMax_=double(max(histImgMax(:)));
            filtMax = double(max(IMGfiltTiled(:)));
            
            if coeffFound % set parameters acc. to the loaded data
                Coeff = CoeffLoaded;
            else
                Coeff = filtMax/5;
            end
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
            isNewTileImage = 1;
            activeSliderNo = 0;
            
            %% LOOP for setting coeff parameters by browsing through frames
            while q == 0 
                if exist('listen1'), delete(listen1);delete(listen2);delete(listen3);delete(listen4); end;
                
                %% update coeff values
                set(hImgCoeffParamText0,'String',Coeff);
                updCoeffFit; % calc coeffFit. (CoeffFit <- hImgCoeffParam)
                if qCoeffSel == 1 % if Coeff selected remove coeff slider
                    delete(hBtnCoeff,hCoeff,hTextCoeff)
                end
                
                figure(hFigCont)
                %% read lookup table settings
                listen1 = addlistener(hCoeff,'ContinuousValueChange',@(hObject, event) updCoeff(hObject, event,hFigDisp,hTextCoeff,hscat,hscat2,hscat3,CoeffFit,sizeImg,nRow,IMGtiled,IMGfiltTiled,filtMax,isDebugFilt,frstFrm,nTile,nFrTile)); 
                listen2 = addlistener(hImgMin,'ContinuousValueChange',@(hObject, event) updImgMin(hObject, event,nTileXY,hTextImgMin,hImg,SHOW,imgMax_, hImgMax, hImgGamma,hLineMin,hLineGamma));
                listen3 = addlistener(hImgMax,'ContinuousValueChange',@(hObject, event) updImgMax(hObject, event,nTileXY, hTextImgMax,hImg,SHOW,imgMax_,hImgMin, hImgGamma,hLineMax,hLineGamma));
                listen4 = addlistener(hImgGamma,'ContinuousValueChange',@(hObject, event) updImgGamma(hObject, event,nTileXY, hTextImgGamma,hImg,SHOW,imgMax_,hImgMin, hImgMax,hLineGamma));
                
                %% browse the stack
                %figure(hFigDisp);
                btn = 0;
                toggleIn = get(hToggleBtnImgWin,'Value'); % Img Win Active
                toggleOut = 0;
                k=0;
                pauseTime = 0;
                while btn == 0 && toggleOut == 0 && toggleIn
                    figure(hFigCont)
                    btn = waitforbuttonpress;
                    if btn == 0, 
                        %pause(2); 
                    end % mouse click, wait for release
                    k = get(hFigCont,'CurrentCharacter');
                    if toggleOut, isNewTileImage = 1; end;
                end
                if toggleOut, k = 0; end
                CoeffSlider = get(hCoeff,'Value'); Coeff = CoeffSlider*filtMax;
                upd = 0; 
                switch lower(k)
                    case 's'
                        CoeffSliderChange =  - 0.001; upd = 1;
                    case 'd'
                        CoeffSliderChange =  + 0.001; upd = 1;
                    case 'w'
                        CoeffSliderChange =  - 0.01; upd = 1;
                    case 'e'
                        CoeffSliderChange =  + 0.01; upd = 1;
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
                    case 'v'
                        nTile = nTile - 100;  isNewTileImage = 1;
                    case 'b' 
                        nTile = nTile + 100; isNewTileImage = 1;
                    case 'n'
                        nTile = nTile - 1500;  isNewTileImage = 1;
                    case 'm' 
                        nTile = nTile + 1500; isNewTileImage = 1;
                end
                if upd, isNewTileImage = 1; end;
                if nTile>floor(nFrames/nFrTile)
                    nTile = floor(nFrames/nFrTile);
                end
                if nTile < 1
                    nTile = 1;
                end             
                
                if upd == 1 % update coeff parameters acc. input keys
                    if activeSliderNo > 0
                        vin = get(hImgCoeffParam(activeSliderNo),'Value'); % value in
                        set(hImgCoeffParam(activeSliderNo),'Value',vin+CoeffSliderChange);
                    else
                        set(hCoeff,'Value',CoeffSlider+CoeffSliderChange); 
                        Coeff = CoeffSlider*filtMax;
                    end
                end       
                if isNewTileImage
                    updCoeffFit; % calc coeffFit. (CoeffFit <- hImgCoeffParam)
                    genPreviewImg; % update image
                    
                    %% DETECTION THRESHOLD
                    CoeffTiled = CoeffFitTiled*Coeff;
                    IMGfilt=IMGfiltTiled; IMG=IMGtiled; CoeffThresh=CoeffTiled;
                    detectThreshold; % (IMGfilt,IMG,CoeffThresh) --> (BW,din)

                    [y,x,v]=find(BW==1);
                    if numel(BW(:)) == sum(BW(:))
                        x=[]; y=[];
                    end;
                    figure(hFigDisp);
                    
                    curCoeff = CoeffFit(framesVec(1))*Coeff;
                    set(gcf,'Name',sprintf('frames:%i - %i, current coeff: %0.2f, suggested coeff: %.02f',1+(nTile-1)*nFrTile, nTile*nFrTile,curCoeff,CoeffSugg))
                    %plot
                    img=SHOW;scaleImg;  % scale image
                    if isDebugFilt, subplot(2,3,1); end
                    hImg = imagesc(img); axis image;
                    hold on;
                    hscat = scatter(x,y,22,'o');hold off;
                    axis image;

                    if isDebugFilt % display debug windows
                        subplot(2,3,2);
                        imagesc(IMGfiltTiled); axis image;
                        hold on;
                        hscat2 = scatter(x,y,22,'o');hold off;

                        subplot(2,3,4);
                        imagesc(IMGfilt1); axis image;
                        hold on;
                        hscat2 = scatter(x,y,22,'o');hold off;

                        subplot(2,3,5);
                        imagesc(IMGfilt2); axis image;
                        hold on;
                        hscat3 = scatter(x,y,22,'o');hold off;
                    IMGfilt2_=IMGfilt2;
                        subplot(2,3,6);
                        imagesc(din); axis image;
                        hold on;
                        hscat3 = scatter(x,y,22,'o');hold off;

                        subplot(2,3,3);
                        imagesc(BW); axis image; 
                    else
                        hscat2 = [];hscat3 = [];
                    end
                    isNewTileImage = 0;                    
                end
                set(hTextCoeff,'String',sprintf('coefficient: %i',round(Coeff)));
                updCoeff(hCoeff, 0,hFigDisp,hTextCoeff,hscat,hscat2,hscat3,CoeffFit,sizeImg,nRow,IMGtiled,IMGfiltTiled,filtMax,isDebugFilt,frstFrm,nTile,nFrTile)
                pause(0.1)
            end % while loop
            CoeffFit = Coeff*CoeffFit; % normalized into CoeffFit
            save('CoeffFit','CoeffFit')
            %save_Coeff;
            close all
            return
            
        end % (if) coeff selection
        close;
        assignFileNames       

        %% find 3D intensity ========================================
        % initialize arrays
        nFrames = Frames-frstFrm2+1;
        time = zeros(nFrames,1); 
        time2 = time;
        ixSptFrm = time;
        errFit7by7 = time;
        nPar = time;
        nSpots = time;
        
        fitVal = zeros(100,nFrames);
        guessNspotsPerFrm=5; % guess for number of spots per frame
        X = zeros(Frames*guessNspotsPerFrm,1);
        X=0;
        Y = X; INT = X; frmNoSpot = X;
        Int = X;
        Px = X; Py = Px;
        
        frameVec = 1:nFrames;

        div = 1;
        
        cfg.fit.Coeff = CoeffFit(1);
        %acqTime = 0.040;
        %numFrmInOneMin = 60/acqTime;
        %numFrmInFiveSec = 60/acqTime;
        i = 1;
        dinPrev2 = 0;dinPrev = 0;
        IMGsum = 0;
        isDoubleHeadedSpot = [];
        %delete('binnedFrames.tif')
        ws = 7; 
        nFramesDisp = 50; 
        IMGmax = 0;
        spotWin = zeros(9);
        h = waitbar(0,'3D localization...'); % ==========================================
        for k=frameVec
            tic
            if exist('BREAK.MAT' ,'file')
                if exist('breakVal.mat','file')
                    clear
                    load breakVal;
                    delete('BREAK.MAT','breakVal')
                else
                    save('breakVal');
                    return;
                end
            end

            frmRead = (k+frstFrm2-1);
            IMG = double(imread(fname,frmRead)); 
            IMG = IMG(yy1:yy2,xx1:xx2,:);
            %imwrite(uint16(IMG),'binnedFrames.tif','WriteMode','append')            
            IMGmean(k) = mean(IMG(:));
            IMGmax(k) = max(IMG(:));

            % DETECTION FILTER
            detectFilter; % IMG --> IMGfilt
            nxt2ndFrm = k+2; if k+2>Frames,nxt2ndFrm=Frames;end;
%CoeffThresh=CoeffFit(1);
            CoeffThresh=CoeffFit(k);
            
            %IMG2 = double(imread(fname,frmRead)); 
            % DETECTION THRESHOLD
            detectThreshold; % (IMGfilt,CoeffThresh) --> (BW,din)
            
            if isempty(find( din > 0, 1)), continue; end;
            %save('filt2','BW','bin','dataFiltDiv','dataFilt','din1','IMG')
            
            [B,L] = bwboundaries(BW,'noholes');
           %% (center of mass localization)                
            nSpots(k) = length(B);
            ixSptFrm(k) = i; % first spot in the frame
            ixSpt = ixSptFrm(k):i-1+nSpots(k); % indices of the spots in the current frame
            for m=1:length(B) % for each spot
                centOfMassLoc; % B --> (X_, Y_)
                if isOptimFit % demagnify
                    X_ = X_/div;
                    Y_ = Y_/div;
                end

                if X_ < 0.5 || X_ > En1-0.5 || Y_ < 0.5 || Y_ > Boy1-0.5 
                    % ignore the spot if in the edge
                    ixSpt(m+1:end) = ixSpt(m+1:end)-1;
                    nSpots(k) = nSpots(k) - 1;
                    continue;
                end

                X(ixSpt(m))=X_;
                Y(ixSpt(m))=Y_;
                INT(ixSpt(m)) = INT_;
                frmNoSpot(ixSpt(m)) = k;

                i = i + 1;
            end
            waitbar(k/Frames,h,'localizing ...')

        end
        ixSptFrm(end+1) = i; % indice of the last spot plus one
        figImgMax = plot(IMGmax);
        print('maxIntensity','-dpng')
        close(h);
        
        close all;

        while ~isempty(find(ixSptFrm==0))
            ixSptFrm(find(ixSptFrm==0))=ixSptFrm(find(ixSptFrm==0)+1);
        end        
        
%        cfg.report.spots_numDetected = numel(xBW);
%        cfg.report.spots_numNeigh = mean(sum(spNonReg7by7>0)+sum(spReg7by7>0));
        disp(sprintf('mean time per frame: %.02f',mean(time)));
        save(xyzDataGausFileNm,'ixSptFrm','X','Y','INT','frmNoSpot','cfg');
        
        %save(spotWinFileNm,'spotWin','xBW','yBW','cfg','spReg7by7','spNonReg7by7','frm'); 
        clear x y x1_ x2_ x1 x2 y1 y2 y1_ y2_ X Y Z INT spotWin NBINs spIn7by7 spReg7by7 spNonReg7by7 spMap frm;
        clear JF J bins JF1 JF2 JF0 yBW xBW Yx Xc BACK
        save(posDataFileNm); % save data 
        
        %% move files
        folderNm= sprintf('%s-coeff%i',label,round(CoeffFit(1)));
        if exist(folderNm,'file')
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
        if exist(coeffMat),copyfile(coeffMat,strcat(folderNm,'/',coeffMat)); end
        if exist('fnameMaxProj'), copyfile(fnameMaxProj,strcat(folderNm,'/',fnameMaxProj));end
        if exist('frames.txt'), copyfile('frames.txt',strcat(folderNm,'/frames.txt'));end
%        copyfile('BWselect.mat',strcat(folderNm,'/','BWselect.mat'));
%        copyfile('CoeffNormMulSmth.mat',strcat(folderNm,'/','CoeffNormMulSmth.mat'));
        %copyfile(fname,strcat(folderNm,'/',fname));
        if exist('spotSel.mat','file'), copyfile('spotSel.mat',strcat(folderNm,'/spotSel.mat')); end;
        if exist('spotSelVal.mat','file'), copyfile('spotSelVal.mat',strcat(folderNm,'/spotSelVal.mat')); end

        movefile(posDataFileNm,strcat(folderNm,'/',posDataFileNm));
        movefile(xyzDataGausFileNm,strcat(folderNm,'/',xyzDataGausFileNm));  
        %copyfile('fname.mat',strcat(folderNm,'/fname.mat'));          
        cd(folderNm)
        save('CFG','cfg');
        fname = ['../' fname ];save('fname','fname');
        delete(traceDataFileNm0)
    end
    %cd ..;    
%return
    
    
    
    
    
    %% TIME To FIND OUT THE TRACES
    if ~exist(traceDataFileNm0,'file')
        load(posDataFileNm);
        if exist(xyzDataGausFileNm,'file')
            load(xyzDataGausFileNm); 
        else
            error('GausFiltFile missing')
            load(xyzDataFileNm); 
        end 
        load CFG; % config
        
        %% fill zero values
        [iy, ix] = find(X>0); %debugSpt = [X(X>0) Y(Y>0) ix];
        while ~isempty(find(ixSptFrm==0))
            ixSptFrm(find(ixSptFrm==0))=ixSptFrm(find(ixSptFrm==0)+1)
        end
        
        %% trace find parameters
        pxSz = str2num(cfg.acq.pixelCAM)/str2num(cfg.acq.magObj)/str2num(cfg.acq.magEx); %[um]
        %sptJmpForTracing = maxTraceSpeed/pxSz*cfg.acq.acqTime;
        
        isFindRecruitment=1;
        sptJmpForTracing = 1; % [px]
        %sptReAppearTime = 1; % no gaps allowed
        sptReAppearTime = 2; 
        minTraceLength = 2;
        traceJmpForCombination = 1;      
                
        cfg.trace = struct;    
        cfg.trace.traceJmpForCombination = traceJmpForCombination;
        cfg.trace.maxTraceSpeed = maxTraceSpeed; 
        cfg.trace.sptJmpForTracing = sptJmpForTracing; % function of maxTraceSpeed
        cfg.trace.sptReAppearTime = sptReAppearTime;
        cfg.trace.minTraceLength = minTraceLength;  
        cfg.trace.minXYspread = minXYspread;
        
        
        %% initialize
        isCropFrames = 0;
        if isCropFrames
            f1=1;f2=300;
            X = X(:,f1:f2);
            Y = Y(:,f1:f2);
        end
        Xilk=X;
        Yilk=Y;
        nSpots = numel(X); % number of spots
        Frames = numel(ixSptFrm)-1; % number of frames
        nTRg = round(nSpots/3); % initial guess for number of traces
        trInf = zeros(nTRg,3); % trace info, memory allocation
        % 1: 1st frame
        % 2: number of frames
        % 3: position in the trace array
        % 4-6: mean x, y , int
        % 7 : std deviation from the center

        %% LOOP
        p=1;
        f=0;
        tic
        isDebug = 0;
        ixTr = 1;
        h = waitbar(0,'Finding the traces...');
        for k=1:Frames-1 % number of frames
            if isDebug, disp(sprintf('=frame#:%i ',k)); end
            Boy = ixSptFrm(k+1) - ixSptFrm(k); % num spots in the frame
            for m_=1:Boy % number of spots
                clear ixSptTr;
%                dif=Inf(Boy,Frames-k+1);
%                difbin=zeros(Boy,Frames-k+1);
                m = m_ + ixSptFrm(k) -1;
                ixSptTr = m;    % first element of the trace
                AslX=X(m);       % last X value in the tracking
                AslY=Y(m);       % last Y value in the tracking
                if AslX == inf
                    continue;
                end
                if isDebug, disp(sprintf('==spot#:%i ',m)); end
                quit = 0;
                ll=k;
                isNewTrace = 0;
                for l=k+1:Frames % later frames
                    if isDebug, disp(sprintf('===check frame#:%i',l)); end
                    if quit, 
                        break; 
                    end
                    ixSptSearchFrm = ixSptFrm(l):ixSptFrm(l+1)-1;
                    dif=sqrt((AslX-X(ixSptSearchFrm)).^2+(AslY-Y(ixSptSearchFrm)).^2);
                    [v n] = min(dif);
                    BOY = Boy;
                    BOY = 1;
                    %for n=1:BOY % all spots
                    if l-ll>sptReAppearTime % if the next frame where a spots re-appears in sptJmp distance is 4 frames apart ignores it.
                        %if sum(sum(difbin(:,l-2:l-1))) == 0, 
                            quit =1;
                            break, 
                        %end
                    end

                    if isDebug, disp(sprintf('===check spot#:%i ',n)); end
                    if dif(n) < sptJmpForTracing; % add to trace
                        %if dif(n,l)==min(dif(:,l))     
                            ll=l;
                            ixSptTr(l-k+1) = ixSptFrm(l)+n-1; % index of the spot in X array
                            AslX=Xilk(ixSptTr(l-k+1));
                            AslY=Yilk(ixSptTr(l-k+1));
                            isNewTrace = 1;
                            if isDebug, disp(sprintf('===add__ spot#:%i, #%i ',n,l-k+1)); end
                            %break
                        %end
                    end
                    %end % spots
                end % frames
                
                if isNewTrace
                    clear tracex tracey traceint;
                    tracex(ixSptTr>0) = Xilk(ixSptTr(ixSptTr>0));
                    tracey(ixSptTr>0) = Yilk(ixSptTr(ixSptTr>0));
                    traceint(ixSptTr>0) = INT(ixSptTr(ixSptTr>0));
                    X(ixSptTr(ixSptTr>0))=Inf;
                    Y(ixSptTr(ixSptTr>0))=Inf;
                    num=numel(find(tracex>0)); % trace length
                    if num>=minTraceLength % if # of data points larger than minTraceLength, than saves as a trace
                        pos=find(tracex>0);
                        ilk=zeros(1,num+1);
                        son=zeros(1,num+1);
                        ilk(1:1:num)=pos(1:1:num);
                        son(2:1:num+1)=pos(1:1:num);
                        fark=ilk-son;
                         %if numel(find(fark==1))>2 % # of consecutive data points
                        p2 = p+numel(tracex)-1; % position of the spot in trace array
                        TraceX(p:p2)= tracex;
                        TraceY(p:p2)= tracey;
                        TraceINT(p:p2)= traceint;   
                        frmNoTrace(p:p2) = k:k+numel(tracex)-1; %% frame numbers of each spot in the trace
                        trInf(ixTr,1) = k; % first frame of the trace
                        trInf(ixTr,2) = numel(tracex); % number of frames (length)
                        trInf(ixTr,3) = p; % position in the trace arrays
                        trInf(ixTr,4) = mean(tracex(tracex>0));
                        trInf(ixTr,5) = mean(tracey(tracex>0));
                        trInf(ixTr,6) = mean(traceint(tracex>0));
                        trInf(ixTr,7) = 0; % std(sqrt((tracex-trInf(ixTr,4)).^2+(tracey-trInf(ixTr,5)).^2)); % RMS deviation
                        trInf(ixTr,8) = findMaxDist(tracex,tracey); % MAX distance
                        %trInf(ixTr,8) = gQF; % gaussian quality factor

                        p=p2+1;
                        ixTr = ixTr + 1; % trace index
                         %end
                    end
                end


            end
           waitbar(k / Frames)
           if exist('_stopRunning-ON','file')
               break
           end
        end
        trInf = trInf(sum(trInf,2)>0,:);
        TrackTime = toc;
%        cfg.report.numSpotInTrace = numel(TraceX(TraceX>0));
%        weight = sum(TraceX>0,2);
%        nHist = hist(weight,max(weight)-min(weight));
%        histTraceLen  = [zeros(min(weight)-1) nHist];
%        cfg.report.histLen = histTraceLen;
        
        save('TrackTime','TrackTime');
        close(h)
        %return; 
        save(traceDataFileNm0,'TraceX','TraceY','TraceINT','trInf','frmNoTrace','cfg'); 
        save(cfgTraceFileNm,'sptJmpForTracing','sptReAppearTime','minTraceLength');
        save('CFG','cfg');
        delete(traceDataFileNm)
    end

    if ~exist(traceDataFileNm) && isCombTraces % speed and trace combination
        %% COMBINE TRACES
        load CFG; % config
        traceJmpForCombination = cfg.trace.traceJmpForCombination;
        load(traceDataFileNm0)
        [m n]=size(TraceX);
        h = waitbar(0,'Combining the traces...');
        disp('Combining the traces...');
        load('xyzDataGaus-coeff941_000.mat', 'ixSptFrm')
        Boy2 = size(trInf,1);
        Frames = numel(ixSptFrm)-1;
        tic;
        com=0;
        i=0;
        while i < size(trInf,1)-1 % each trace
            i=i+1;
        %     TraceDif1=0;TraceDif2=0;TraceDif3=0;
            LastElement = trInf(i,3)+trInf(i,2)-1; % (position in the trace array)+(number of frames)
            LastBefore=LastElement-1;
            j=i;
            while j < size(trInf,1)
                j=j+1;
                FirstElement=trInf(j,3); % (position in the trace array)
                FirstAfter=FirstElement+1;
                diffTime = FirstElement-LastElement; % [frames]
                if (0 <= diffTime)  && (diffTime <= 2) 
                    TraceDif1=sqrt([TraceX(LastElement)-TraceX(FirstElement)]^2+[TraceY(LastElement)-TraceY(FirstElement)]^2);
                    TraceDif2=sqrt([2*TraceX(LastElement)-TraceX(LastBefore)-TraceX(FirstElement)]^2+[2*TraceY(LastElement)-TraceY(LastBefore)-TraceY(FirstElement)]^2);
                    TraceDif3=sqrt([TraceX(LastElement)-2*TraceX(FirstElement)+TraceX(FirstAfter)]^2+[TraceY(LastElement)-2*TraceY(FirstElement)+TraceY(FirstAfter)]^2);
                    TraceDif=[TraceDif1,TraceDif2,TraceDif3];
                    %TraceDif=sqrt((TraceX(i,LastElement)-TraceX(j,FirstElement))^2+(TraceY(i,LastElement)-TraceY(j,FirstElement))^2);
                    if min(TraceDif)<traceJmpForCombination % combine traces
                        isOverlapping = 0;
                        isSpace = 0;
                        if diffTime == 2 % overlapping
                            isSpace = 1; % one frame is empty 
                        elseif diffTime == 0 % overlapping
                            isOverlapping = 1;
                            k1 = LastElement;
                            k2 = FirstElement;
                            TraceX(k1)= [TraceX(k1)*TraceINT(k1)+TraceX(k2)*TraceINT(k2)]/[TraceINT(k1)+TraceINT(k2)];
                            TraceY(k1)= [TraceY(k1)*TraceINT(k1)+TraceY(k2)*TraceINT(k2)]/[TraceINT(k1)+TraceINT(k2)];
                            TraceZ(k1)= [TraceZ(k1)*TraceINT(k1)+TraceZ(k2)*TraceINT(k2)]/[TraceINT(k1)+TraceINT(k2)];
                            TraceINT(k1)= [TraceINT(k1)*TraceINT(k1)+TraceINT(k2)*TraceINT(k2)]/[TraceINT(k1)+TraceINT(k2)];
                        end
                        trInf(i,2) = trInf(i,2)+trInf(j,2)-isOverlapping+isSpace; % number of frames
                        % join the traces
                        fr1   = LastElement+1+isSpace; % extension of the first trace
                        fr2   = LastElement+isSpace+trInf(i,2);
                        fr1_2 = FirstElement+isOverlapping; % coord. of the second trace
                        fr2_2 = FirstElement+trInf(j,2)-1;
                        TraceX(fr1:end) = [TraceX(fr1_2:fr2_2) TraceX(fr1:fr1_2-1) TraceX(fr2_2+1:end) ];

                        %frmNoTrace(fr1:end) = 
                        LastElement=trInf(i,1)+trInf(i,2)-1;
                        LastBefore=LastElement-1;
                        trInf(j:end-1,:) = trInf(j+1:end,:);
                        trInf = trInf(1:end-1,:);

                        % update trace info of the following traces
                        trInf(j:end,3)=trInf(j:end,3)-isOverlapping+isSpace;
                        trInf(i+1:j-1,3)=trInf(i+1:j-1,3) + fr2_2-fr1_2+isSpace;
                        % update trace values

                        %trInf(i,4) = mean(tracex(tracex>0));
                        %trInf(i,5) = mean(tracey(tracex>0));
                        %trInf(i,6) = mean(traceint(tracex>0));
                        %trInf(i,7) = std(sqrt((tracex-trInf(ixTr,4)).^2+(tracey-trInf(ixTr,5)).^2)); % RMS deviation                            

                        % update 

                        com=com+1;
                    end
                end
            end % i <= Boy2-1
        end

        save('CFG','cfg');
        save(traceDataFileNm,'TraceX','TraceY','TraceZ','TraceINT','TraceSpeed','trInf','frmNoTrace')
        delete(traceJmplessDataFileNm)
    end
    

    %return
    %% PLOT1 : fill jumps
    if ~exist(traceJmplessDataFileNm,'file')
        load CFG; % config
        disp('preparing plot data...')
        load(traceDataFileNm0);
        
        load(xyzDataGausFileNm, 'ixSptFrm')
        if exist(xyzDataGausFiltFileNm)
            load(xyzDataGausFiltFileNm);
        elseif exist(xyzDataGausFileNm)
            load(xyzDataGausFileNm); 
        else
            load(xyzDataFileNm); 
        end
        % trInf :
        % 1: first frame of the trace
        % 2: number of frames (length)
        % 3: position in the trace arrays
        
        
        tic
        hWBfillgaps = waitbar(0,'filling gaps');
        TraceX2 = TraceX;
        TraceY2 = TraceY;        
        for tr = 1:size(trInf,1) % number of traces
            %% fill gaps (missing frames in jumpy traces)
            
            traceX = TraceX(trInf(tr,3):trInf(tr,3)+trInf(tr,2)-1);
            traceY = TraceY(trInf(tr,3):trInf(tr,3)+trInf(tr,2)-1);
            traceX(traceX<0) = 0.001; % no negative values
            traceY(traceY<0) = 0.001;
            
            ind = find(traceX~=0);
            frst = min(ind); last = max(ind);
            %iJump = find(trace(frst:last)==0)+frst-1;

            jmp = abs((traceX(frst:last)>0)-1);
            if sum(jmp > 0)
                bnd = bwboundaries(jmp,'noholes');
                for i = 1:numel(bnd) % for each jump
                    temp = bnd{i}+frst-1; % jump boundaries
                    jb = temp(:,2); clear temp;
                    mx= max(jb); mn=min(jb);
                    jL = mx-mn+2; % length

                    jSx = traceX(mx+1)-traceX(mn-1);% size
                    jSy = traceY(mx+1)-traceY(mn-1);% size
                    jsX = jSx/jL;% step
                    jsY = jSy/jL;% step
                    jVx = traceX(mn-1)+jsX*(1:jL-1);
                    jVy = traceY(mn-1)+jsY*(1:jL-1);
                    traceX(mn:mx)=jVx;
                    traceY(mn:mx)=jVy;
                    TraceX2(trInf(tr,3):trInf(tr,3)+trInf(tr,2)-1) = traceX;
                    TraceY2(trInf(tr,3):trInf(tr,3)+trInf(tr,2)-1) = traceY;
                end
            end  
            waitbar(tr/size(TraceX,1))
        end
        close(hWBfillgaps)
        toc
        save(traceJmplessDataFileNm,'TraceX2','TraceY2')
    end

%return % un comment to stop





    
    
    
    
    % w = warning('query','last');id = w.identifier;warning('off',id)
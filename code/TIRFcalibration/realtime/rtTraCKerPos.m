function rtTraCKerPos
% RUN in waSeq\tracker\
% 'style','text','BackgroundColor',[1 1 1],
    % F = findall(0,'type','figure','tag','TMWWaitbar'); delete(F)
pause(2)
cd('E:\MATLAB\TIRFcalibration\data\Ata01_5_125X100Y50x50_realtime')    
    cd waSeq\tracker\
    cfg = '..\..\cfgRT';
    cfg_=load(cfg);
    cfg = cfg_.cfg;
    ndigit = cfg.ndigit; % # of digits for sequence number
    En1 = cfg.w;
    Boy1 = cfg.h;
    label = cfg.label;
          
    tic; logFN = cfg.logPos; fid = fopen(logFN,'w'); wait = 0;
    clck = clock; fprintf(fid,'start time m= %2i secs=%6.03f\n',clck(5),clck(6));
    
    %% detection parameters
    WindowSize = 3; 
    BigWindowSize=WindowSize+2;
    
    lap=[-1,-1,-1;-1,8,-1;-1,-1,-1];
            gausKernelSz = 3;
            gausKernelSg = 0.7;
            gausKernelSzPadded = 5;
    pdSz = (gausKernelSzPadded - gausKernelSz)/2; % pad size
    gausHat1=fspecial('gaussian', gausKernelSz, gausKernelSg);

    elev=1/(gausKernelSzPadded^2-gausKernelSz^2);        
    if elev==inf, elev=0; end
    gausHat1 = gausHat1 + elev;
    gausHatPad1 = padarray(gausHat1,[pdSz pdSz]);
    gausHatPad1 = gausHatPad1 - elev;
    gausHat1 = gausHatPad1;        

        
	
    %% read filename       
    fn_ = load('..\..\fname0.mat'); % filename_WA_
    fname = fn_.fname0;
    digitFormat = sprintf('''%%0%1ii''',ndigit);
    

    %% feedback to calling function
    fcall = 'rtTraCKerPos';
    fdbck.inWait = 0;
    fdbck.inWaitCounting = 0;
    fdbck.inPause = 0;
    fdbck.inSave = 0;
    fdbck.inSaveCounting = 0;
    fdbck.inSaveCountingIX = 0;
    fdbck.inSaveCountingMAX = cfg.inSaveCountingMAX;
    fdbck.inStop = 0;    
    IMG = zeros(50);
    IMGfilt = IMG;
    BW = IMG;
    din = IMG;
    n = 1;
    while (1)
        time = toc; fprintf(fid,'while loop n=%3i time=%6.03f\n',n,time);
        
        while (1) % wait for update
            coeffFN = dir('Coeff*.mat');
            if isempty(coeffFN), pause(0.01); continue; end
            Coeff_=load(coeffFN.name);
            Coeff = Coeff_.Coeff;
            if numel(Coeff) < n
                if rem(wait,10) == 0, time = toc; fprintf(fid,'wait for   n=%3i time=%6.03f, ncoeff:%i %s\n',n,time,numel(Coeff),coeffFN.name); wait = 1; end
                wait = wait + 1;
                %[fdbck] = funcFeedback(cfg.msgTXT,fdbck,fcall);
                %if fdbck.inStop, break;  end % STOP
            else
                time = toc; wait = 0; fprintf(fid,'updated    n=%3i time=%6.03f\n',n,time);
                break; % continue
            end 
            pause(0.010)
        end
        try
            Coeff = Coeff(n);
            digitTXT = eval(['sprintf(' digitFormat ',n)'] );
            fnameSeq = ['..\' fname 'WA_' digitTXT '.tif'];
            CoeffThresh = Coeff(end);
            posDataFileNm = sprintf('posData-coeff%i_%s_%s.mat',round(Coeff),label,digitTXT);

            %% FIND X & Y
            IMG = double(imread(fnameSeq,1)); % READ IMAGE
    %            IMG(90:94,1:4)=IMG(95:99,1:4); % ======remove
            IMGfilt = detectFilter; % IMG --> IMGfilt
            [din,BW] = detectThreshold; % (IMGfilt,CoeffThresh) --> (BW,din)
            if isempty(find( din > 0, 1)), continue; end
            [B,~] = bwboundaries(BW,'noholes');
            dbg = 0;
            if dbg
                dbgDetectFN = 'dbgDetect.tif';
                delete(dbgDetectFN);
                imwrite(uint16(IMG),dbgDetectFN,'tif','WriteMode','append')
                imwrite(uint16(IMGfilt),dbgDetectFN,'tif','WriteMode','append')
                imwrite(uint16(din),dbgDetectFN,'tif','WriteMode','append')
                imwrite(uint16(BW),dbgDetectFN,'tif','WriteMode','append')
            end

            %% (center of mass localization)
            ixSpt = 1;
            X = []; Y = []; INT = [];
            for m=1:length(B) % for each spot
                [X_,Y_] = centOfMassLoc; % B --> (X_, Y_)

                if X_ < 0.5 || X_ > En1-0.5 || Y_ < 0.5 || Y_ > Boy1-0.5 
                    % ignore the spot if in the edge
                    continue;
                end
                X(ixSpt)=X_;
                Y(ixSpt)=Y_;
                INT(ixSpt) = INT_;
                ixSpt = ixSpt + 1;
            end 
            save(posDataFileNm,'X','Y','INT');
            n = n + 1;
        catch ME
            
            time = toc; fprintf(fid,'ERROR: %s\n',ME.identifier);
        end
            
        
    end
    
    
    
%===============================================================
    function IMGfilt = detectFilter
        
        IMGfilt2 = nan(size(IMG));
        IMGfilt1 = imfilter(IMG,gausHat1,'symmetric');
        IMGfilt1 = imfilter(IMGfilt1,lap,'symmetric');
        IMGfilt = max([IMGfilt1(:) IMGfilt2(:)],[],2);
        IMGfilt = reshape(IMGfilt,size(IMGfilt1));
    end
%===============================================================
    function [din,BW] = detectThreshold
        % thresholding
        % (IMGfilt,IMG,CoeffThresh) --> (BW,din)
        dataFiltDiv = IMGfilt./CoeffThresh;
        bin = im2bw(dataFiltDiv,1);
        din = uint16(bin).*uint16(IMG);
        BW = imregionalmax(din, 8);
    end

%===============================================================
    function [X_,Y_] = centOfMassLoc
        k=1;
        % center of mass localization
        % B --> (X_, Y_, INT_)
        % B: boundaries of the BW peaks
        % (X_, Y_): coordinates of the localized peak
        c=cell2mat(B(m));
        Py=uint16(mean(c(:,1)));
        Px=uint16(mean(c(:,2)));
        Px0 = double(Px);
        Py0 = double(Py);

        % adjust the big window center position for the spots at
        % the edges
        PxBW = Px;
        PyBW = Py;

        if (Px-(BigWindowSize+1)/2)<1
            PxBW=(BigWindowSize+1)/2;
        end
        if (Py-(BigWindowSize+1)/2)<1
            PyBW=(BigWindowSize+1)/2;
        end
        if (Px+(BigWindowSize+1)/2)>En1
            PxBW=En1-(BigWindowSize+1)/2;
        end
        if (Py+(BigWindowSize+1)/2)>Boy1
            PyBW=Boy1-(BigWindowSize+1)/2;
        end
        dPy = double(PyBW)-double(Py);
        dPx = double(PxBW)-double(Px);
        Px = PxBW; Py = PyBW;
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
        BACK =min(BACKmean);

        %FIND Total Intensity
        INT_=sum(sum(Window))-BACK * (WindowSize)^2;


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

        for ii=1:WindowSize
           TopY(ii)=sum(Window(ii,:));
        end
        TopY=TopY-min(TopY);
        TopColum=sum(TopY);

        for ii=1:WindowSize
            WSumY=WSumY+ii*TopY(ii);
        end

        Xc(k)=WSumX/TopRow;
        Yc(k)=WSumY/TopColum;
        if isnan(Xc(k)), Xc(k)=(WindowSize+1)/2; end
        if isnan(Yc(k)), Yc(k)=(WindowSize+1)/2; end

        PXc=uint8(Xc(k));
        PYc=uint8(Yc(k));

        %center of intensity
        X_=double(Px)+Xc(k)-double((WindowSize+1)/2);
        Y_=double(Py)+Yc(k)-double((WindowSize+1)/2);

        X_ = X_-dPx;
        Y_ = Y_-dPy;
    end
   
end
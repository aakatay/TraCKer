function rtTraCKerPos(varargin)
% RUN in waSeq\tracker\
% 'style','text','BackgroundColor',[1 1 1],
    % F = findall(0,'type','figure','tag','TMWWaitbar'); delete(F)
    
    isCallOutside = 0;
    if nargin == 1
        isCallOutside = 1;
        cfg = varargin{1};
        PWD = cfg.pwd;
        cd(PWD)
        q=0;
        pause(1);
        while q == 0
            try 
                cd waSeq\tracker\
                q = 1;
            catch
                q=0;
                %save PWD PWD
                %pause(1);
            end
        end
    else
        cfg = '..\..\cfgRT';    cfg_=load(cfg);    cfg = cfg_.cfgSave;
    end

    waWin = cfg.waWin; % walkiong average window length
    ndigit = cfg.ndigit; % # of digits for sequence number
    En1 = cfg.w;
    Boy1 = cfg.h;
    label = cfg.label;
    tloopPause = cfg.tloopPause;
    isTlog = cfg.isTlog;      
    timeOut = cfg.timeOut;    
    tic; 
    if isTlog, logFN = cfg.logPos; fid = fopen(logFN,'w'); wait = 0; end
    if isTlog, clck = clock; fprintf(fid,'start time m= %2i secs=%6.03f\n',clck(5),clck(6)); end
    
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
    btnMAT                  = '..\..\signals\btnMAT.mat';
    MATrtTraCKerPos         = '..\..\signals\MATrtTraCKerPos.mat';
    quitToutMAT             = '..\..\signals\quitTout.mat';
    syncFrameMAT            = '..\..\signals\syncFrame.mat';
    % fdbck inputs to funcFeedback : 
    fdbck.nFrst         = 0;
    fdbck.nLast         = 0;
    fdbck.runProcess    = 0;
    fdbck.syncWait      = 0;
    fdbck.toutOn        = 0;
    fdbck.syncHere      = 0;
    fdbck.isStop        = 0;
    fdbck.ssSnap        = 0;
    fdbck.ssSave        = 0;
    % fdbck inputs/outputs to funcFeedback: 
    fdbck.isSS          = 0;
    fdbck.inSS          = 0;
    fdbck.dispSS        = 0;
    
    IMG = zeros(50);
    IMGfilt = IMG;
    BW = IMG;
    din = IMG;
    n = 1;
    isStop = 0;
    tout =[]; % timeout
    syncFrameList = [];
    while (1)
        if isTlog, time = toc; fprintf(fid,'while loop n=%3i time=%6.03f\n',n,time); end
        
        while (1) % wait for update
            coeffFN = dir('Coeff*.mat');
            if isempty(coeffFN), pause(0.01); continue; end
            Coeff_=loadMAT(coeffFN.name); Coeff = Coeff_.Coeff;

            b_=tryLoadbtnMAT(sprintf('out=load(''%s'');',btnMAT),cfg.tTryLoop); btnStart = b_.btnStart; btnSync = b_.btnSync; btnSnap = b_.btnSnap; btnSave = b_.btnSave; btnStop = b_.btnStop; 
            if btnStop >= 0
                isStop = 1;
                break;
            end
            if fdbck.toutOn == 1 % timeout
                if exist(quitToutMAT) % quit timeout
                    fdbck.toutOn = 0;
                else
                    continue;
                    if isTlog, if wait == 0, time = toc; fprintf(fid,'timeout@   n=%3i time=%6.03f\n',n,time);wait = 1;end; end
                end
            elseif fdbck.toutOn == -1
                if fdbck.runProcess % reset timeout
                    fdbck.toutOn = 0;
                    tout = [];
                end   
            end
            
            if fdbck.syncWait % wait for sync
                while ~exist(syncFrameMAT) % wait for sync data
                    pause(0.01)
                end
                if fdbck.syncHere
                    if 1 || btnSync == -1 % reset sync
                        fdbck.syncWait = 0;
                    end
                else
                    syncMAT=load(syncFrameMAT); nLastSync = syncMAT.nLast; % sync frame
                    n = nLastSync;
                    syncFrameList = [syncFrameList n];
                    fdbck.syncHere=1; 
                end
            elseif btnSync >= 0 && btnStart>=1 % check new sync
                if exist(syncFrameMAT)
                    syncMAT=load(syncFrameMAT); nLastSync = syncMAT.nLast; % sync frame
                    if ~ismember(nLastSync,syncFrameList) % new sync
                        fdbck.syncWait = 1;
                    end
                else % new sync
                    fdbck.syncWait = 1;
                end 
            end            
            fdbck.nFrst = n;
            fdbck.nLast = n + waWin - 1;
            
            fdbck.runProcess = 0;
            if numel(Coeff) >= n % newData
                if fdbck.syncHere, fdbck.syncHere = 0; fdbck.syncWait = 0; end
                fdbck.runProcess = 1; % process update
                tout = toc; % reset timeout time
            elseif fdbck.toutOn==0
                if isempty(tout)
                    tout = toc; % time wait
                elseif toc-tout > timeOut % timeout 
                    fdbck.toutOn = 1;
                end
            end    
            
            [fdbck] = funcFeedback(cfg,fdbck,fcall);
            if fdbck.runProcess % process new data
                if isTlog, wait = 0; time = toc; fprintf(fid,'updated    n=%3i time=%6.03f\n',n,time); end
                clck = clock;                 
                break; 
            else % wait
                if isTlog, if wait == 0, time = toc; fprintf(fid,'wait for   n=%3i time=%6.03f\n',n,time);wait = 1;end; end
            end 
            pause(tloopPause)
        end
        if isStop
            if exist('fid'), fclose(fid);end
            lmpState = -1; save(MATrtTraCKerPos,'lmpState','-append'); % stop
            break;
        end
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
        
    end
    
    
    
%===============================================================

    function CoeffOut = loadMAT(coeffFN)
        ME = 1;
        while ~isempty(ME)
            ME = [];
            try 
                CoeffOut=load(coeffFN);
            catch ME
                pause(0.001)
            end
        end
    end
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
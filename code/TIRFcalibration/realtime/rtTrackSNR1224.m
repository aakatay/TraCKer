% reads the 'filename_001.tif' and data from tracking of 'filename_002.tif' 
% processes 'filename_001.tif'
function rtTrackSNR(varargin)
%cd('E:\MATLAB\TIRFcalibration\data\Ata01_5_125X100Y50x50_realtime');    
    isdbgAcquisition = 1; % generated input files for debug
    
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
                cd waSeq\tracker\rtData
                q = 1;
            catch
                q=0;
                %save PWD PWD
                %pause(1);
            end
        end
    else
        %% load cfgRT
        cfg = '..\..\..\cfgRT';
        c_=load(cfg);
        cfg = c_.cfgSave;
    end
    
    dbgSel = 0;
    dbgSNRvoronIMG = 0;
    dbgSNRimg = 0; % SNR movie
    dbgSNR = 0; % intensity traces
    SNRmul = 1000; % for intesity scale 'SNRimg.tif'
    nt0 = 3000; % # of traces to initialize
   
    
    waWin       = cfg.waWin; % walkiong average window length
    ndigit      = cfg.ndigit; % # of digits for sequence number
    label       = cfg.label;
    w           = cfg.w;
    h           = cfg.h;
    wsz         = cfg.wsz; % window size for SM crop
    szXY        = [w h];
    szYX        = fliplr(szXY);
    acqTime     = cfg.acqTime; % [s]
    stdWin      = cfg.stdWin; % number of frames to calc. std
    mag         = cfg.dispMag; % display size
    scrnSzIn    = cfg.scrnSzIn;
    mag         =cfg.dispMag;
    SNRcolorThresh = cfg.SNRcolorThresh;
    tloopPause  = cfg.tloopPause;
    isTlog      = cfg.isTlog;
    timeOut     = cfg.timeOut;
    tic;
    if isTlog, logFN = cfg.logSNR; fid = fopen(logFN,'w'); c = onCleanup(@()fclose(fid)); wait = 0; end
    if isTlog, clck = clock; fprintf(fid,'start time m= %2i secs=%6.03f\n',clck(5),clck(6)); end

    digitFormat = sprintf('%%0%1ii',ndigit);
    s = floor(wsz/2); 
    scrnSzIn = [1600 900];
    scrnSz =get(0,'ScreenSize');
    scrnSz = scrnSz(3:4);
    scrnSzLeftBottom = scrnSz-scrnSzIn;
    [mag, pos, szx, szy ] = calcMaxMag(zeros(szXY),mag,scrnSzIn);
    szXYmag = [szx, szy];
    szYXmag = fliplr(szXYmag);
    SNRmovVoronoi = zeros([szYXmag 3]);
    
    %% input files
    traceFN = 'traceData_';
    fn_ = load('..\..\..\fname0'); % filename_WA_
    fname = fn_.fname0;
    fname = ['..\..\..\acq\' fname];   
    
    if isdbgAcquisition
        inputDIR = '..\..\..\input\';
        inputFN_ = dir([inputDIR '*.tif']);
        inputFN = ['..\..\..\input\' inputFN_(1).name];
        inputIx = 2;
        clckAcqTime0 = 0;
    end
        
    if isCallOutside
        %% Communication variables
        % MAT files assoc. to parallel workers
        MATrtWAmean         = '..\..\..\signals\MATrtWAmean.mat';
        MATrtDetectThresh   = '..\..\..\signals\MATrtDetectThresh.mat';
        MATrtTraCKerPos     = '..\..\..\signals\MATrtTraCKerPos.mat';
        MATrtTraCKerTrace   = '..\..\..\signals\MATrtTraCKerTrace.mat';
        MATrtTrackSNR       = '..\..\..\signals\MATrtTrackSNR.mat';
        MATcell = {MATrtWAmean,MATrtDetectThresh,MATrtTraCKerPos,MATrtTraCKerTrace,MATrtTrackSNR};
        lmp1 = cfg.lmps{1}; % rtWAmean
        lmp2 = cfg.lmps{2}; % rtDetectThresh
        lmp3 = cfg.lmps{3};
        lmp4 = cfg.lmps{4}; % rtTraCKerTrace
        lmp5 = cfg.lmps{5};


        btnStart0 = cfg.btns.btnStart0;
        btnSync0 = cfg.btns.btnSync0;
        btnSnap0 = cfg.btns.btnSnap0;
        btnSave0 = cfg.btns.btnSave0;
        btnStop0 = cfg.btns.btnStop0;
        BTNcell0 = {btnStart0,btnSync0,btnSync0,btnSave0,btnStop0};
        
        
        btnCol = [cfg.btnColDefault;cfg.btnColPress;cfg.btnColActive];
        
        isQuit = [];
        i = [];
    end
    %% output files
    SNRmovieFN          = '..\..\..\snapSaveOUT\SNRmovie.tif';
    SNRmovieConvFN      = '..\..\..\snapSaveOUT\SNRmovieConv.tif';
    SNRmovieVoronoiFN0  = '..\..\..\snapSaveOUT\SNRmovieVoronoi'; % .tif label
    SNRmovieVoronoiFN   = [];
    lensConfigFN0       = '..\..\..\snapSaveOUT\lensConfig'; % .txt label
    SNRplotFN           = '..\..\..\snapSaveOUT\SNRplot.tif';
    numSMplotFN         = '..\..\..\snapSaveOUT\SNRnumSMplot.tif';
    snrSTDplotFN        = '..\..\..\snapSaveOUT\SNRstdplot.tif';
    %SNRdataFN = 'SNRdata.mat';
    
    %% SNR image figure
    if dbgSNRvoronIMG
        tit = 'SNR voronoi image';
        %pos(1) = pos(1) - 800;
        %pos(1) = pos(1)- 700;        
        figSNRvoronIMG = figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos szXYmag(1) szXYmag(2)]);
        axeSNRvoron = axes('Parent',figSNRvoronIMG,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren','XLim',0.5+[0 szXYmag(1)],'YLim',0.5+[0 szXYmag(2)]);
    end
    
    %% SNR plot
    tit = 'SNR plot';
    szx2 = szx;
    szy2 = scrnSzIn(2)-szy-45;
    pos2 = [pos(1) scrnSzLeftBottom(2)];
    figSNRplot = figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos2 szx2 szy2]);
    %axeSNRplot = axes('Parent',figSNRplot,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren','XLim',0.5+[0 szXYmag(1)],'YLim',0.5+[0 szXYmag(2)]);

    if dbgSNRimg
        tit = 'SNR image';
        pos(1) = pos(1) - 800;
        pos(1) = pos(1)- 700;            
        figSNRimg = figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos/2 szXYmag(1) szXYmag(2)]);
        axeSNR = axes('Parent',figSNRimg,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren','XLim',0.5+[0 szXYmag(1)],'YLim',0.5+[0 szXYmag(2)]);
    end

    %% crop
    xc1 = 160;
    yc1 = szXY(2)-110;
    xc2 = xc1+20;
    yc2 = yc1+25;

    xc1 = 1;
    yc1 = szXY(2)-1;
    xc2 = xc1+100;
    yc2 = yc1+100;

    xc1 = 1;
    yc1 = 1;
    xc2 = xc1+50;
    yc2 = yc1+50;
    
    %% feedback to calling function
    fcall = 'rtTrackSNR';
    btnMAT          = '..\..\..\signals\btnMAT.mat';
    MATrtTraCKSNR   = '..\..\..\signals\MATrtTrackSNR.mat';
    quitToutMAT     = '..\..\..\signals\quitTout.mat';
    
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
    % others
    
    
    %% display figure
    if dbgSNR
        figSNR = figure(1212);
        dispSz = [1400 900];
        tit = 'Trace Plots';
        set(figSNR,'DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),...
            'Position',[250 50 dispSz(1) dispSz(2)]);
    end
    
    %% first frame
    fnameSeq = [fname num2str(1,digitFormat) '.tif'];
    Afrst = imread(fnameSeq);
    Afrst = padarray(Afrst,[s s]);
    
    %% loop
    n = 1;
    f = 1; % a index
    aIX = cell(nt0,1);
    fIX = cell(nt0,1);
    mlast = 0; % SM index for SNR
    ixFrm = 0; % frst frame (0 traces)
    XYS = [];
    snrMean = [];
    tout =[]; % timeout
    while (1)        
        if isTlog, time = toc; fprintf(fid,'while loop n=%3i time=%6.03f\n',n,time); end
        timeLoop(n) = toc;
        %% load data
        while (1) % wait for update
            twhile1 = toc;
            traceSeq = [traceFN label '_' num2str(n+1,digitFormat) '.mat'];
            
            b_=load(btnMAT); btnStart = b_.btnStart; btnSync = b_.btnSync; btnSnap = b_.btnSnap; btnSave = b_.btnSave; btnStop = b_.btnStop; 
            if btnStop >= 0
                fdbck.isStop = 1;
                if isQuit 
                    break;
                else % wait updateCommunication 
                    continue; 
                end
            end
            if fdbck.toutOn == 1 % timeout
                if exist(quitToutMAT) % quit timeout
                    fdbck.toutOn = -1;
                else
                    continue;
                    if isTlog, if wait == 0, time = toc; fprintf(fid,'timeout@   n=%3i time=%6.03f\n',n,time);wait = 1;end; end
                end
            elseif fdbck.toutOn == -1
                if fdbck.runProcess % reset timeout
                    fdbck.toutOn = 0; tout = [];
                    delete(quitToutMAT)
                end                
            end
            
            if fdbck.syncWait % wait for sync
                while ~exist(syncFrameMAT) % wait for sync data
                    pause(0.01)
                end
                syncMAT=load(syncFrameMAT); nLastSync = syncMAT.nLast; % sync frame
                n = nLastSync;
                if fdbck.syncHere % reset sync
                    fdbck.syncWait = 0; 
                    fdbck.syncHere = 0;
                    delete(syncFrameMAT);
                    % push btnSync=-1
                else
                    fdbck.syncHere=1; 
                end
            elseif btnSync >= 0 && btnStart==1 % enter sync
                fdbck.syncWait = 1; % push btnSync=1
            end    
            
            if btnSnap == 0 || btnSave == 0 % snap/save
                if (ssSnap || ssSave) && dispSS % reset snap/save
                    ssSnap = 0; ssSave = 0;
                elseif isSS == 0 && fdbck.syncHere == 1 % enter Snap/Save
                    isSS = 1;
                    if btnSnap == 0,fdbck.ssSnap = 1;elseif btnSave == 0,fdbck.ssSave = 1;end                
                    setSNRmovieVoronoiFN;
                end
            end

            fdbck.nFrst = n+1;
            fdbck.nLast = n + waWin;
            
            fdbck.runProcess = 0;
            if exist(traceSeq) % newData
                if fdbck.syncWait 
                    if fdbck.syncHere % process update
                        fdbck.runProcess = 1;
                    end
                else % process update
                    fdbck.runProcess = 1; 
                end
                tout = toc; % reset timeout time
            elseif fdbck.toutOn==0 % default
                if isempty(tout) % run timer
                    tout = toc; % time wait
                elseif toc-tout > timeOut*10 % timeout 
                    fdbck.toutOn = 1;
                end
            end    
            
            [fdbck] = funcFeedback(cfg,fdbck,fcall);
            updateCommunication;
            dbgAcquisition; % emulate camera acquisition
            if fdbck.runProcess % process new data
                if isTlog, wait = 0; time = toc; fprintf(fid,'updated    n=%3i time=%6.03f\n',n,time); end
                clck = clock;                 
                break; 
            else % wait
                if isTlog, if wait == 0, time = toc; fprintf(fid,'wait for   n=%3i time=%6.03f\n',n,time);wait = 1;end; end
                %writeSNRvoronoiMov
                %playSNRvoronoiMov
                %runAnalysisWin
                if 0
                    figure;
                    timeLoop = timeLoop(2:end)-timeLoop(1:end-1);
                    tSM = tSM(2:end);
                    plot([timeLoop;tSM]');
                    title(sprintf('mean time : %.03f',mean(timeLoop(9:end))));
                    legend('loop' ,'SM')
                end
            end 
            twhile2 = toc; twhileDiff = tloopPause - (twhile2 - twhile1); pause(twhileDiff)
        end      
        if fdbck.isStop
            if exist('fid'), fclose(fid);end
        end
        n = n + 1;
        
        if 0 && n == 2 % pauses the rtWAmean to sync
            test = []; save test test
            fopen('..\..\..\logData\syncSignal.txt','w');
        end

        q=0;
        while q == 0
            try 
                fn_ = load(traceSeq); % trace MAT
                TraceX = fn_.TraceX;
                TraceY = fn_.TraceY;
                trInf = fn_.trInf;
                q = 1;
            catch
                q=0;
                %save PWD PWD
                pause(.01);
            end
        end        

        %% read trace data
        X = TraceX;
        Y = TraceY;
        
        % SM selection
        %trInf = trInf(515:518,:);
        %frm1 = trInf(:,1);
        %trInf(frm1~=1)=[]; % remove new recruitments

        frm1 = trInf(:,1);
        frm2 = trInf(:,2)+trInf(:,1)-1;
        nt = size(trInf,1); % # of single molecule
        if size(aIX)<nt, aIX{nt} = []; end % increase array size
        if size(fIX)<nt, fIX{nt} = []; end

        dbg = 0;
        if dbg    
            dbgImg = zeros(szXYmag);
            xdbg = round(trInf(:,4));
            ydbg = round(trInf(:,5));
            dbgImg(sub2ind(size(dbgImg),ydbg,xdbg)) = 1;
            figure(209); imagesc(dbgImg); axis image
            imwrite(dbgImg,'dbgImg.tif')
        end        

        %% PROCESS FRAME ===========================================================
        if dbgSNRimg, snrIMG = zeros(szXYmag); end
        
        fnameSeq = [fname num2str(n,digitFormat) '.tif'];
        A = imread(fnameSeq);
        A = padarray(A,[s s]);
        IX = find((frm1<=n) .* (n<=frm2));
        fTrace = n - frm1(IX) + 1;
        
        x = ceil(cellfun(@(v) v(end), X(IX)));
        y = ceil(cellfun(@(v) v(end), Y(IX)));    
        
        %% select isolated SM
        t1 = toc;
        findSMisolated; % (x,y) -> ixs
        t2 = toc;
        tSM(n) = t2-t1;
        ixSM = IX(ixs);
        if isempty(ixSM), continue; end

        %% crop single molecule windows
        
        tt =0;
        % construct crop image stack
        for j = 1:numel(ixSM) % each single molecule in that frame
            ixsm = ixSM(j); % single molecule index
            y_ = y(ixs(j));
            x_ = x(ixs(j));
            %if ((x_<xc1) || (y_<yc1) || (x_>xc2) || (y_>yc2)), tt = tt + 1;continue; end
            if ((x_<=s) || (y_<=s) || (x_+s-1>szXY(1)) || (y_+s-1>szXY(2))), continue; end % at the edge
            if n == 2 % add the first frame data
                a(:,:,f) = Afrst(y_-s:y_+s,x_-s:x_+s); % crop
                aIX{ixsm} = [aIX{ixsm} f]; % index for the frame of the single molecule
                fIX{ixsm} = [fIX{ixsm} 1]; % frames
                f=f+1;
            end
            a(:,:,f) = A(y_-s:y_+s,x_-s:x_+s); % crop
            aIX{ixsm} = [aIX{ixsm} f]; % index for the frame of the single molecule
            fIX{ixsm} = [fIX{ixsm} n]; % frames
            f=f+1;
        end
        
        %% calculate SNR
        m = 0;
        Xs = []; Ys = []; Frm = []; TRinf = []; B = []; INT = []; SNR = []; MIX = [];
        if n >= stdWin 
            aIXix_ = find(~cellfun(@isempty,aIX)); % all SM found so far
            aIXix = aIXix_;
            aIXix( (   (n-frm1(aIXix)+1)<stdWin   )) = []; % remove SM shorter than stdWin
            aIXix( (   frm2(aIXix) < n  )) = []; % remove SM bleach before the current frame (n)
            if dbgSNRimg, snrImg = zeros(szXYmag); end
            num_ixImg2(n) = 0;
            for j = 1:numel(aIXix) % for each SM 
                ix = aIXix(j);
                fr = fIX{ix}; % frames
                if fr(end) ~= n
                    if dbgSel, disp('1nodata'); end
                    continue; end % no data in this frame
                if numel(fr)-stdWin+1 < 1,if dbgSel,disp('2notenough');end; continue; end % not enough data
                if fr(end)-stdWin+1 ~= fr(end-stdWin+1) % missing data in the last stdWin frames
                    if dbgSel,disp('3missingdata'); end
                    continue;
                end
                m = m + 1;
                num_ixImg2(n) = num_ixImg2(n) + 1;
                aix = aIX{ix};
                asm = a(:,:,aix); % SM crop stack
                smLast = asm(:,:,end); % last image

                intWin = double(smLast(s-1:s+2,s-1:s+2)); % 4by4
                int = sum(intWin(:));
                peak = max(intWin(:));
                asm(s-1:s+2,s-1:s+2,:) = nan;
                a2D = double(reshape(asm,wsz^2,size(asm,3)));

                astd = std(a2D,0,1);
                astd(isnan(astd)) = [];
                b = mean(astd); % background std
                snr = peak/sqrt(peak+b^2);

                Xs(m) = ceil(cellfun(@(v) v(end), X(ix)));
                Ys(m) = ceil(cellfun(@(v) v(end), Y(ix)));
                Frm(m) = n;
                TRinf(m,:) = trInf(ix,:);     
                B(m) = b;
                INT(m) = int;
                SNR(m) = snr;                   
                MIX(m) = ix; % SM index
                
                
                if dbgSNRimg
                    Xmag = round(Xs(m)*mag);
                    Ymag = round(Ys(m)*mag);
                    snrImg(Ymag,Xmag,j) = SNR(m);
                end

                %BCK(n) = mean(bckgrnd(:));
                %INT(n) = max(peak(:));
                %RECT(n,:) = rect;
            end
            
            if dbgSNRimg, snrIMG = sum(snrImg,3); end
            %sn{n} = find(snrImg~=0);
            %% display
            if dbgSNR
                figure(figSNR)
                runIntensityTracesRT
            end
            %% save
            %save(SMdataFN,'Xs','Ys','Frm','TRinf','B','INT','SNR','MIX','aIX','fIX');
        end
        
        if dbgSNRimg
            figure(figSNRimg);
            set(axeSNR,'Parent',figSNRimg,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','XLim',1+[0 szXYmag(1)],'YLim',1+[0 szXYmag(2)]);
            imagesc(flipud(snrIMG))
            pause(acqTime)
            SNRimg(:,:,n) = snrIMG;
        end
        
        %sum(reshape(SNRimg,260^2,226),1)

        if ~isempty(Xs)
            XYS = [XYS;Xs' Ys' SNR']; 
        end
        %save(SNRdataFN,'ixFrm','XYS')
        
        
        %% display
        mlast = mlast+m;
        ixFrm = [ixFrm mlast]; % indices of spots in each frame
        
        % SNR data
        isContinue = 0;
        xys = XYS(ixFrm(n-1)+1:ixFrm(n),:);
        if isempty(xys)
            SNr =[]; isContinue = 1;
            snrMean(n) = nan;
        else
            SNr = xys(:,3);
            snrMean(n) = max(SNr);
        end
        
        % SNR voron
        if fdbck.inSS || dbgSNRvoronIMG
            figure(figSNRvoronIMG)
            t1 = toc;
            SNRmovVoronoi = writeSNRvoronoiFrame;
            t2 = toc;
            %figure(figSNRvoronIMG)
            %imagesc(SNRmovVoronoi);
            tVoron(n) = t2-t1;
            %pause
        end
        
        if fdbck.dispSS % show mean snap/save image
            
            
        end
        
        if isContinue, continue; end
        
        % SNR plot
        figure(figSNRplot)
        nBars = 13;
        snrBarsDisp = snrMean;
        xd1 = 1;
        if n > nBars
            xd1 = n - nBars+1;
            %snrBarsDisp = snrMean(xd1:n);
        end
        bar(snrBarsDisp)
        L1 = SNRcolorThresh(1);
        L2 = SNRcolorThresh(2);
        xd2 = xd1 + nBars;
        line([0 xd2],[L1 L1],'Color','r')
        line([0 xd2],[L2 L2],'Color','g')
        xlim([xd1-0.5 xd2])
        ylim([0 10])
        xticks(xd1:xd1+nBars-1);
        

        if fdbck.isStop, break; end
        
        continue
         
%% STOP

        if dbg
            ix = find(~cellfun(@isempty,aIX));
            ix(1:numel(MIX),2)=MIX';

        end

        if 0  % write SM stack
            %%
            ix = find(~cellfun(@isempty,aIX));
            ixs=ix(1);
            frm = aIX{ixs};
            a(:,:,frm);
            xs = round(trInf(ixs,4));
            ys = round(trInf(ixs,5));
            fname = sprintf('SM%03i_x%iy%i',ixs,xs,ys);
            stackWrite(a,fname)

        end
    end % while
    
    
    
%% FUNCTIONS --------------------------------------------------------    

    %% COMMUNICATION vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    function dbgAcquisition
        % emulate camera acquisition
        if ~isdbgAcquisition, return; end
        try
            IMGin = imread(inputFN,inputIx);
        catch
            isdbgAcquisition = 0; return;
        end
        inputFN2  = sprintf('..\\..\\..\\acq\\%s_%04i.tif',inputFN(16:end-4),inputIx);
        clckAcq = clock; clckAcqTime = clckAcq(4)*24+clckAcq(5)*60+clckAcq(6);
        if clckAcqTime - clckAcqTime0 > cfg.acqTime
            IMGin = imread(inputFN,inputIx);
            imwrite(IMGin,inputFN2);
            inputIx = inputIx + 1;
        end
    end

    function updateCommunication
        if ~isCallOutside, return; end
        isQuit = 0;
            %btnSave = 1; save(btnMAT,'btnSave','-Append');
            %btnSnap = 1; save(btnMAT,'btnSnap','-Append');

        if ~isequal(get(btnStart0,'BackgroundColor'), cfg.btnColPress), return;  end % wait for startBtn
        checkStop; % if btnStop
        if isQuit, return; end
        for i = 1:5 % read inputs from five parallel workers
            matFN = MATcell{i};
            matfn = tryLoad(sprintf('out=load(''%s'');',matFN),cfg.tTryLoop); % lmpState n

            % update lamp displays
            lmp = cfg.lmps{i};
            if matfn.lmpState==-1 % stop
                set(lmp,'Color',cfg.lmpColStop); 
            elseif matfn.lmpState==-1i % syncWait
                set(lmp,'Color',cfg.lmpColSyncWait); 
                if i == 5, set(btnSync0,'BackgroundColor',cfg.btnColActive); end % send syncReady 
            elseif matfn.lmpState==0 % timeout
                set(lmp,'Color',cfg.lmpColTimeout); 
            elseif matfn.lmpState==1i % syncHere
                set(lmp,'Color',cfg.lmpColSyncHere); 
            elseif matfn.lmpState==1 % active
                checkResetSync;
                set(lmp,'Color',cfg.lmpColActive); 
            end

            % update frame numbers
            set(cfg.frmnoFrst{i},'Text',sprintf('%04i',matfn.nFrst)); % update text labels
            set(cfg.frmnoLast{i},'Text',sprintf('%04i',matfn.nLast));
        end

        % update button states
        BTNcell = {btnStart,btnSync,btnSnap,btnSave,btnStop}; % variables
        if btnStart == 0 && n==1, btnStart=1; end % ready for start sync
        for j = 1:5 % read inputs from five button press
            set(BTNcell0{j},'BackgroundColor', btnCol(BTNcell{j}+2,:));
        end
        
        
        pause(0.01)

    end

    function checkResetSync
        if i ~= 5, return; end
        if isequal(get(lmp,'Color'), cfg.lmpColSyncHere) % reset sync
            set(btnSync0,'BackgroundColor',cfg.btnColDefault);
        end
        
    end


    function checkStop
        % check if btnStop
        if ~isequal(get(btnStop0,'BackgroundColor'), cfg.btnColPress), return; end
        % check if all workers stopped
        if ~isequal(get(lmp1,'Color'), cfg.lmpColStop), return; end
        if ~isequal(get(lmp2,'Color'), cfg.lmpColStop), return; end
        if ~isequal(get(lmp3,'Color'), cfg.lmpColStop), return; end
        if ~isequal(get(lmp4,'Color'), cfg.lmpColStop), return; end
        if ~isequal(get(lmp5,'Color'), cfg.lmpColStop), return; end
        % => all stopped
        set(btnStop0,'BackgroundColor',cfg.btnColActive);
        btnStop = 1; save(btnMAT,'btnStop','-Append');
        isQuit = 1;
    end
    %% COMMUNICATION -----------------------------------------------------------









    function findSMisolated 
        % find non-overlapping SM
        smMap = zeros(szYX);
        smMap(sub2ind(szYX,y,x)) = 1;
        smMapCv = conv2(smMap,ones(wsz),'same');
        smMapCv(smMapCv~=1) = -10000;
        smMapCv2 = conv2(smMapCv,ones(wsz),'same');
        smMapCv3 = zeros(szYX);
        smMapCv3(smMapCv2 == wsz^2) = 1;
        %smMapCv4 = circshift(smMapCv3,[1 1]);
        smMapCv4 = smMapCv3;
        smMapCv5 = smMapCv4+smMap;
        ixImg = find(smMapCv5==2);
        ixIMG = sub2ind(szYX,y,x)';
        ixs = find(ismember(ixIMG,ixImg));
        

        dbg2 = 0;
        if dbg2 
            figure(344);imagesc(smMap)
            figure(345);imagesc(smMapCv)
            figure(346);imagesc(smMapCv2)
            figure(347);imagesc(smMapCv3)
            %figure(348);imagesc(smMapCv4)
            figure(349);imagesc(smMapCv5)
            cc=4;
        end
    end
        
    function runSNRmov
    end

    function runSNRplot
        %% SNR plot
        figSNRplot = figure(201);
        SNR2d = reshape(SNRimg,[szXYmag(2)*szXYmag(1) nfr]);
        SNR2dSUM = sum(SNR2d,1);
        SNR2dNUM = sum(SNR2d>0,1);
        snrmean = SNR2dSUM./SNR2dNUM;

        plot(snrmean); 
        title('SNR mean')
        xlabel('frames')
        ylabel('SNR')   
        imgFig = getframe(gcf);
        imgOut = imgFig.cdata;
        imwrite(imgOut,SNRplotFN);
    end
    
    function runSNRstdPlot
        %% SNR std plot
        figSNRstdPlot = figure(202);
        SNR2d = reshape(SNRimg,[szXYmag(2)*szXYmag(1) nfr]);
        SNR2dSTD = std(SNR2d,1);

        plot(SNR2dSTD); 
        title('SNR STD')
        xlabel('frames')
        ylabel('SNR STD')   
        imgFig = getframe(gcf);
        imgOut = imgFig.cdata;
        imwrite(imgOut,snrSTDplotFN);
    end

    function runNumSMplot
        %% number of SM plot
        figNumSMplot = figure(204);
        plot(num_ixImg2);
        title('number of SM for SNR reading')
        xlabel('frames')
        ylabel('Number of Single Molecules')   
        imgFig = getframe(gcf);
        imgOut = imgFig.cdata;
        imwrite(imgOut,numSMplotFN);
    end    

    function runIntensityTracesRT
        %% intensity traces
        hold on;
        %load(SMdataFN); % 'Xs','Ys','Frm','TRinf','B','INT','SNR','MIX');
        smix = unique(MIX); % # selected SM
        if isempty(smix), return; end
        if smix(1)==0, smix(1)=[]; end % remove zero
        ns = numel(smix);

        intShft = 10e4;
        intShft = 0;
        ixs = [];
        for i = 1:ns
            ix = smix(i); % SM index
            mix = find(MIX==ix);
            frm = Frm(mix);
            FRM(i)=numel(frm);
            %if numel(frm)<217,continue; end
            ixs = [ixs ix];
            int = INT(mix);
            INTsave(:,i) = int;
            plot(frm,int+intShft*(i-1));
        end
        save('SNRintTraces','INTsave','frm','ixs')
        savefig('SNRintTraces');
        hold off;
    end

    function runIntensityTraces
        %% intensity traces
        figIntensityTraces = figure(1212); maximize;clf;
        hold on;
        %load(SMdataFN); % 'Xs','Ys','Frm','TRinf','B','INT','SNR','MIX');
        smix = unique(MIX); % # selected SM
        if smix(1)==0, smix(1)=[]; end % remove zero

        smixSel = [11 14 16 21];
        ns = numel(smix);
        ns2 = 0;
        ns3 = 0;
        FRM = [];
        intShft = 10e4;
        intShft = 0;
        ixs = [];
        for i = 1:ns
            ixc = smix(i); % SM index
            mix = find(MIX==ixc);
            frm = Frm(mix);
            FRM(i)=numel(frm);
            if numel(frm)<100,continue; end
            if numel(frm)<217,continue; end

            ns2 = ns2 + 1;
            %if ~ismember(ns2,smixSel),continue; end;
            ns3= ns3+1;
            ixs = [ixs ixc];
            int = INT(mix);
            INTsave(:,ns3) = int;
            subplot(4,1,ns3)
            plot(frm,int+intShft*(ns3-1));
            ylim([0 max(int)*1.2])
            pause(0.5)
        end
        save('SNRintTraces','INTsave','frm','ixs')
        savefig('SNRintTraces');
        hold off;
    end


%% ============ Analysis ===================================================
    function runAnalysisWin
        fig= uifigure('Position',[20 200 300 300]);

        x1 = 20;
        dy = 30;
        y1 = [0:3]*dy+10;

        w = 100;
        h1 = 20;
        btn1 = uibutton(fig,'Position',[x1 y1(1) w h1],'Text','Intensity','ButtonPushedFcn', @(btn1,event) dispTraceInt);
        btn2 = uibutton(fig,'Position',[x1 y1(2) w h1],'Text','SNR','ButtonPushedFcn', @(btn2,event) runPause);
        btn3 = uibutton(fig,'Position',[x1 y1(3) w h1],'Text','SNR movie','ButtonPushedFcn', @(btn3,event) writeSNRmov);
        btn4 = uibutton(fig,'Position',[x1 y1(4) w h1],'Text','SNR voronoi mov','ButtonPushedFcn', @(btn4,event) writeSNRvoronoiMov);
        while(1)
            pause(1)
        end
    end

    function dispTraceInt
        figure(901)
        hold on
        smix = unique(MIX); % # selected SM
        if isempty(smix), return; end
        if smix(1)==0, smix(1)=[]; end % remove zero
        ns = numel(smix);

        intShft = 10e4;
        intShft = 0;
        ixs = [];
        for i = 1:ns % each selected trace
            ix = smix(i); % SM index
            mix = find(MIX==ix);
            frm = Frm(mix);
            FRM(i)=numel(frm);
            %if numel(frm)<217,continue; end
            ixs = [ixs ix];
            int = INT(mix);
            %INTsave(:,i) = int;
            plot(frm,int+intShft*(i-1));
        end
        save('SNRintTraces','INTsave','frm','ixs')
        savefig('SNRintTraces');      
        hold off  
        
    end


    function writeSNRvoronoiMov
        delete(SNRmovieVoronoiFN);
        
        CM = jet(256);
        CM = gray(256);
        for i = 2:n % frame index
            writeSNRvoronoiFrame(i);
        end
    end


    function SNRmovVoronoi = writeSNRvoronoiFrame
        CM = gray(256);
        if numel(SNr) < 5 % write a blank image
            SNRmovVoronoi = zeros(szXY*mag);
            %imwrite(SNRmovVoronoi,SNRmovieVoronoiFN,'WriteMode','append','Compression', 'none') 
        else
%SNRcolorThresh = [.1 5];
            % SNR intensity range scales
            if max(SNr)<SNRcolorThresh(1)
                sscale = SNRcolorThresh(1);
                EdgeColorSel = 1; % red
            elseif max(SNr)<SNRcolorThresh(2)
                sscale = SNRcolorThresh(2);
                EdgeColorSel = 2; % green 
            elseif max(SNr)<SNRcolorThresh(3)
                sscale = SNRcolorThresh(3);
                EdgeColorSel = 3; % blue
            else % max > SNRcolorThresh(3)
                sscale = SNRcolorThresh(3)*2;
                EdgeColorSel = 3; % blue
                
            end
            [SNRmovVoronoi,tPatch(n)] = getVoronoinImg(figSNRvoronIMG,xys(:,1:2),SNr/sscale,szXY,mag,CM,EdgeColorSel);
            imwrite(SNRmovVoronoi,SNRmovieVoronoiFN,'WriteMode','append','Compression', 'none') 
        end
    end

    function setSNRmovieVoronoiFN
        if btnSnap % change filename for multiple snaps
            SNRmovieVoronoiFNheading = ['snap' SNRmovieVoronoiFN0];
            lensConfigFNheading = ['snap' lensConfigFN0];
            fdir=rdir([SNRmovieVoronoiFNheading '*.tif']); 
            [~,flast]=max(cell2mat({fdir.datenum})); 
            ixFN = 1; % first snap
            if ~isempty(flast) % another file
                flastFN = fdir(flast).name;
                ixFN = str2num(flastFN(end-5:end-4))+1;
            end
            SNRmovieVoronoiFN = [SNRmovieVoronoiFNheading sprintf('%02i.tif',ixFN)];
            lensConfigFN = [lensConfigFNheading sprintf('%02i.txt',ixFN)];
        elseif btnSave 
            SNRmovieVoronoiFN = ['save' SNRmovieVoronoiFN0 '.tif'];
            lensConfigFN = ['save' lensConfigFN0 '.tif'];
        end
        
        % save lens config
        lensParam = {'L1tilt' 'L1shft' 'L1dist' 'L2tilt' 'L2shft' 'L2dist'};
        if feedbck.inSS==1 % first snap
            fid = fopen(lensConfigFN,'w');
            fprintf(fid,'first frame:%i\n',n);
            for i=1:6 % each paramter
                lensCfg(i) = get(cfg.lensBox(i),'Value');
                fprintf(fid,'%s:%.02f\n',lensParam{i},lensCfg(i));
            end
            fclose(fid);
        end
    end
        

    function playSNRvoronoiMov
        fnMov = dir([SNRmovieConvFN(1:end-4) '*.tif']);
        if ~isempty(fnMov)
        end
    end

    function writeSNRmov
        if ~dbgSNRimg, return; end
        %% SNR movie
        if ~isempty(find(SNRimg*SNRmul>=2^16)), warning('SNR image saturated use lower a SNRmul');end
        %SNRmovieConvFN2 = [ SNRmovieConvFN(1:end-4) sprintf('frm%04i-%04i.tif',) ]
        SNRmov = SNRimg*SNRmul;
        stackWrite(SNRmov,SNRmovieFN);
        cvWin = ones(5);
        SNRmovCv = convn(SNRmov,cvWin,'same');
        stackWrite(SNRmovCv,SNRmovieConvFN);
    end

end
































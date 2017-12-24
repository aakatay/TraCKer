function rtWAmean(varargin)
% walking average mean
% RUN in mainfolder
    if nargin == 1
        cfg = varargin{1};
    else
        cfg_=load('cfgRT'); cfg = cfg_.cfgSave;
    end
    
    waWin = cfg.waWin; % walkiong average window length
    acqTime = cfg.acqTime; % [s]
    ndigit = cfg.ndigit; % # of digits for sequence number
    isTlog = cfg.isTlog;
    tloopPause = cfg.tloopPause;
    timeOut = cfg.timeOut;
    
    tic;
    if isTlog, logFN = cfg.logWA; fid = fopen(logFN,'w'); wait = 0; end
    if isTlog, clck = clock; fprintf(fid,'start time m= %2i secs=%6.03f\n',clck(5),clck(6)); end
    clck = clock;
        
    digitFormat = sprintf('%%0%1ii',ndigit);
    outDIR = 'waSeq\';

    %% set filenames and array
    fname0_ = load('fname0');
    fname0 = fname0_.fname0;
        % reads the position of crop from the filename
        
    Awa = zeros(cfg.h,cfg.w,waWin);

        cx0 = cfg.crop(1)+1;
        cy0 = cfg.crop(2)+1;
        cx1 = cx0 + cfg.w -1;
        cy1 = cy0 + cfg.h -1;

    %% walking average
    n = 1;
    
    %% feedback to calling function
    fcall = 'rtWAmean';
    btnMAT              = 'signals\btnMAT.mat';
    MATrtWAmean         = 'signals\MATrtWAmean.mat';
    quitToutMAT         = 'signals\quitTout.mat';
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
    
    
    nSync = waWin+2;
    isStop = 0;
    tout =[]; % timeout
    while (1)
        %pause(0.2)
        if 0 && n == nSync % wait for sync signal
            while ~exist('syncSignal.txt')
                D = dir([fname0 '*.tif']);
                CD = struct2cell(D);
                dates = cell2mat(CD(6,:));
                [~,ix]=sort(dates);
                lastFn = D(ix(end)).name;
                lastn = str2num(lastFn(end-ndigit-3:end-4));
                nJump = lastn - nSync;
                save('logData\syncJump','nJump');
                n = lastn;
            end
        end
        
        if isTlog, time = toc; fprintf(fid,'while loop n=%3i time=%6.03f\n',n,time); end
        %tl(n) = toc;
        while (1) % wait for update
            dbgWaitAcqTime(clck)
            fnameSeq = ['acq\' fname0 num2str(n,digitFormat) '.tif'];
            b_=load(btnMAT); btnStart = b_.btnStart; btnSync = b_.btnSync; btnSnap = b_.btnSnap; btnSave = b_.btnSave; btnStop = b_.btnStop; 
            if btnStop >= 0
                isStop = 1;
                break;
            end
            if fdbck.toutOn == 1 % timeout
                if exist(quitToutMAT) % quit timeout
                    fdbck.toutOn = 0;
                    save(quitToutMAT); % send quitTimeout
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
                if fdbck.syncHere
                    if btnSync == -1 % reset sync
                        fdbck.syncWait = 0;
                        fdbck.syncHere = 0;
                    end
                elseif btnSync==1
                    fdir=rdir([fname0 '*.tif']); [~,flast]=max(cell2mat({fdir.datenum})); fnameSeq = fdir(flast).name; 
                    n=sscanf(fnameSeq,[fname0 '%f.tif']);
                    fdbck.syncHere=1; save(syncFrameMAT,'nLast'); % process new image
                end
            elseif btnSync >= 0 && btnStart==1
                fdbck.syncWait = 1;
            end
            fdbck.nFrst = n - waWin + 1;
            fdbck.nLast = n;
            
            fdbck.runProcess = 0;
            if exist(fnameSeq) % newData
                if fdbck.syncWait 
                    if fdbck.syncHere % process update
                        fdbck.runProcess = 1;
                    end
                else % process update
                    fdbck.runProcess = 1; 
                end
                tout = toc; % reset timeout time
            elseif fdbck.toutOn==0
                if isempty(tout)
                    tout = toc; % time wait
                elseif toc-tout > timeOut % timeout 
                    fdbck.toutOn = 1;
                end
if n>waWin, pause(2); end
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
            lmpState = -1; save(MATrtWAmean,'lmpState','-append'); % stop
            break;
        end
        
        Awa = circshift(Awa,[0 0 1]);
        Awa_ = imread(fnameSeq);
        Awa(:,:,1) = Awa_(cy0:cy1,cx0:cx1);
        
        
        if n>=waWin
            A = uint16(mean(Awa,3));
            nWA = n - waWin+1;
            fnameWaSeq = [outDIR fname0 'WA_' num2str(nWA,digitFormat) '.tif'];
            ME = 1;
            while ~isempty(ME)
                ME = [];
                try imwrite(A,fnameWaSeq);
                catch ME
                    disp(sprintf('error in writing, trying again | frame number: %i',i)); %#ok<DSPS>
                end
                
                %pause(acqTime/10)
            end
        end
        n = n + 1; % # of frames
        %pause(acqTime/3)
        
        inWait = 0;
    end
    n = n - 1;
    
    function dbgWaitAcqTime(clck)
        %pauses for acqTime
        %return
        sec0 = clck(5)*60+clck(6);
        while (1)
            clck = clock;
            sec = clck(5)*60+clck(6);
            if sec-sec0>=acqTime
                break
            end
        end
    end












    function timelapsePlot
        return;
        figure;
        tloop = tl(2:end)-tl(1:end-1);
        plot([ tloop;twrt ]')
        legend('loop' ,'write')
    end
end

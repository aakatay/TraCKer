function rtDetectThresh(varargin)
%% CORRECTION: remove a line cropping an area
% detects threshold for tracking
% output Coeff.mat
% RUN in waSeq\


% image =======
% filtered image =======
% filtered image peak MASK 
% image PEAKs =======
% histogram peak (image PEAKs)

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
                cd waSeq\
                q = 1;
            catch
                q=0;
                %save PWD PWD
                %pause(1);
            end
        end
    else
        load ..\cfgRT;
        cfg = cfgSave;
    end        
    
    
    %test = []; save test test
    
    dbg = 0;

    
    waWin = cfg.waWin; % walkiong average window length
    acqTime = cfg.acqTime; % [s]
    ndigit = cfg.ndigit; % # of digits for sequence number
    w = cfg.w;
    h = cfg.h;
    bc = cfg.bc; % [nM] background concentration
    tloopPause = cfg.tloopPause;
    isTlog = cfg.isTlog;
    timeOut = cfg.timeOut;
    tic; 
    if isTlog, logFN = cfg.logThresh; fid = fopen(logFN,'w'); wait = 0; end
    if isTlog, clck = clock; fprintf(fid,'start time m= %2i secs=%6.03f\n',clck(5),clck(6)); end
    

    digitFormat = sprintf('%%0%1ii',ndigit);
    outDIR = 'tracker\';
    
    CoeffFN = [outDIR 'Coeff_' cfg.label '.mat'];
    Coeff = []; save(CoeffFN,'Coeff');

    % coefficient multiplier
    cm1 = bc/10;
    cm2 = bc/27;
    coeffMul1 = 1+0.4; 
    coeffMul2 = 1-0.15; 
    coeffMul1 = 1+cm1; % corrects shift of histogram peak to the right
    coeffMul2 = 1-cm2; % corrects for increased intensity fluctuations

    coeffMul = coeffMul1*coeffMul2; % 4nM background


    %% set filenames and array
    load('..\fname0.mat'); % % filename_WA_
    
    %% feedback to calling function
    fcall = 'rtDetectThresh';
    btnMAT              = '..\signals\btnMAT.mat';
    MATrtDetectThresh   = '..\signals\MATrtDetectThresh.mat';
    quitToutMAT         = '..\signals\quitTout.mat';
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

    n = 1;
    isStop = 0;
    tout =[]; % timeout
    while (1)
        if isTlog, time = toc; fprintf(fid,'while loop n=%3i time=%6.03f\n',n,time); end

        while (1) % wait for update
            fnameSeq = [fname0 'WA_' num2str(n,digitFormat) '.tif'];

            b_=load(btnMAT); btnStart = b_.btnStart; btnSync = b_.btnSync; btnSnap = b_.btnSnap; btnSave = b_.btnSave; btnStop = b_.btnStop; 
            if btnStop >= 0
                isStop = 1;
                break;
            end
            if fdbck.toutOn == 1 % timeout
                if exist(quitToutMAT) % quit timeout
                    fdbck.toutOn = 0;
                else
                    continue;
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
                syncMAT=load(syncFrameMAT); nLastSync = syncMAT.nLast; % sync frame
                n = nLastSync;
                if fdbck.syncHere
                    if btnSync == -1 % reset sync
                        fdbck.syncWait = 0;
                        fdbck.syncHere = 0;
                    end
                else
                    fdbck.syncHere=1; 
                end
            elseif btnSync >= 0 && btnStart==1
                fdbck.syncWait = 1;
            end           
            fdbck.nFrst = n;
            fdbck.nLast = n + waWin - 1;
            
            fdbck.runProcess = 0;
            if exist(fnameSeq) % newData
                if fdbck.syncWait 
                    if fdbck.syncHere % process update
                        fdbck.runProcess = 1;
                    end
                else % process update
                    fdbck.runProcess = 1; 
                end
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
            lmpState = -1; save(MATrtDetectThresh,'lmpState','-append'); % stop
            break;
        end

        %% image
        IMG = imread(fnameSeq,1);
%                    IMG(90:94,1:4)=IMG(95:99,1:4); % ======remove
        [~,c] = loghist(IMG,30,0);
        hstep = c(2)-c(1); % hist step

        if dbg
            nc = 6; % # columns
            figure(505);
            subplot(2,nc,1)
            imagesc(invert(IMG));
            colormap('gray');axis image; title('image')

            subplot(2,6,nc+1)
            [~,c] = loghist(IMG,30,1);
        end

        %% filtered image 
        gausKernelSz = 3;
        gausKernelSg = 0.7;
        gausKernelSzPadded = 5;

        pdSz = (gausKernelSzPadded - gausKernelSz)/2; % pad size
        gausHat1=fspecial('gaussian', gausKernelSz, gausKernelSg);
        lap=[-1,-1,-1;-1,8,-1;-1,-1,-1];

        IMGfilt2 = nan(size(IMG));
        IMGfilt1 = imfilter(IMG,gausHat1,'symmetric');
        IMGfilt1 = imfilter(IMGfilt1,lap,'symmetric');
        IMGfilt = max([IMGfilt1(:) IMGfilt2(:)],[],2);
        IMGfilt = reshape(IMGfilt,size(IMGfilt1));
        if dbg
            subplot(2,nc,2)
            imagesc(invert(IMGfilt))
            colormap('gray');axis image;title('filtered image')
            imwrite(IMGfilt,'IMGfilt.tif')
            subplot(2,6,nc+2)
        end

        [~,c] = loghist(IMGfilt,30,dbg);

        %% filtered image peak MASK 
        BW = imregionalmax(IMGfilt, 8);

        %% image PEAKs
        IMGpeak = uint16(BW).*IMG;
        if dbg
            subplot(2,nc,3)
            imagesc(invert(IMGpeak))
            colormap('gray');axis image; title('image PEAKs')
            imwrite(IMGpeak,'IMGpeak.tif')

            subplot(2,6,nc+3)
            [hc,c] = loghist(IMGpeak,50,dbg); % hist count and centers
            xlabel('intensity'); ylabel('counts')
        end
        %[hc,c] = loghist(IMGpeak,50,1); % hist count and centers


        %%  histogram peak(image PEAKs)
        [hc,c] = loghist(IMGpeak,50,0); % hist count and centers
        [~,locs] = findpeaks(hc);
        peakInt1 = c(locs(1));

        Coeff(n) = peakInt1*coeffMul;
        save(CoeffFN, 'Coeff');
        n = n + 1; % # of frames
        
    end
    n = n - 1;
end
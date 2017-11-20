function rtWAmean
% walking average mean
% RUN in mainfolder

    
    cfg_=load('cfgRT');
    cfg = cfg_.cfg;
    
    
    waWin = cfg.waWin; % walkiong average window length
    acqTime = cfg.acqTime; % [s]
    ndigit = cfg.ndigit; % # of digits for sequence number
    digitFormat = sprintf('%%0%1ii',ndigit);
    outDIR = 'waSeq\';
    if exist(outDIR), rmdir(outDIR,'s'); end
    mkdir(outDIR)

    %% set filenames and array
    fname0_ = load('fname0');
    fname0 = fname0_.fname0;
        % reads the position of crop from the filename
        
    Awa = zeros(cfg.h,cfg.w,waWin);


    %% walking average
    n = 1;
    
    %% feedback to calling function
    fcall = 'rtWAmean';
    fdbck.inWait = 0;
    fdbck.inWaitCounting = 0;
    fdbck.inPause = 0;
    fdbck.inSave = 0;
    fdbck.inSaveCounting = 0;
    fdbck.inSaveCountingIX = 0;
    fdbck.inSaveCountingMAX = cfg.inSaveCountingMAX;
    fdbck.inStop = 0;
    
    while (1)
        while (1) % wait for update
            fnameSeq = [fname0 num2str(n,digitFormat) '.tif'];
            if ~exist(fnameSeq) % wait for update
                [fdbck] = funcFeedback(cfg.msgTXT,fdbck,fcall);
                if fdbck.inStop, break;  end % STOP
            else
                break; % continue
            end 
        end
        
        Awa = circshift(Awa,[0 0 1]);
        Awa_ = imread(fnameSeq);
        cx0 = cfg.crop(1)+1;
        cy0 = cfg.crop(2)+1;
        cx1 = cx0 + cfg.w -1;
        cy1 = cy0 + cfg.h -1;
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
                pause(acqTime/10)
            end
        end
        n = n + 1; % # of frames
        pause(acqTime/3)
        
        inWait = 0;
    end
    n = n - 1;
    
end

function rtWAmean
% walking average mean
% RUN in mainfolder
    cfg_=load('cfgRT');
    cfg = cfg_.cfg;
        
    tic; logFN = cfg.logWA; fid = fopen(logFN,'w'); wait = 0;
    clck = clock; fprintf(fid,'start time m= %2i secs=%6.03f\n',clck(5),clck(6));
    
    waWin = cfg.waWin; % walkiong average window length
    acqTime = cfg.acqTime; % [s]
    ndigit = cfg.ndigit; % # of digits for sequence number
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
    fdbck.inWait = 0;
    fdbck.inWaitCounting = 0;
    fdbck.inPause = 0;
    fdbck.inSave = 0;
    fdbck.inSaveCounting = 0;
    fdbck.inSaveCountingIX = 0;
    fdbck.inSaveCountingMAX = cfg.inSaveCountingMAX;
    fdbck.inStop = 0;
    
    nSync = waWin+2;
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
        
        time = toc; fprintf(fid,'while loop n=%3i time=%6.03f\n',n,time);
        tl(n) = toc;
        while (1) % wait for update
            dbgWaitAcqTime(clck)
            fnameSeq = [fname0 num2str(n,digitFormat) '.tif'];
            if ~exist(fnameSeq) % wait for update
                if 0
                    figure;
                    tloop = tl(2:end)-tl(1:end-1);
                    plot([ tloop;twrt ]')
                    legend('loop' ,'write')
                end
                if wait == 0, time = toc; fprintf(fid,'wait for   n=%3i time=%6.03f\n',n,time);wait = 1;end
                [fdbck] = funcFeedback(cfg.msgTXT,fdbck,fcall);
                if fdbck.inStop, break;  end % STOP
            else
                wait = 0; time = toc; fprintf(fid,'updated    n=%3i time=%6.03f\n',n,time);
                clck = clock;
                break; % continue
            end 
            %pause(0.1)
        end
        
        
        Awa = circshift(Awa,[0 0 1]);
        Awa_ = imread(fnameSeq);
        Awa(:,:,1) = Awa_(cy0:cy1,cx0:cx1);
        
        
        if n>=waWin
                t1 = toc;
            A = uint16(mean(Awa,3));
                t2 = toc;
                twrt(n) = t2-t1;
            nWA = n - waWin+1;
            fnameWaSeq = [outDIR fname0 'WA_' num2str(nWA,digitFormat) '.tif'];
            ME = 1;
            while ~isempty(ME)
                ME = [];
                t1 = toc;
                try imwrite(A,fnameWaSeq);
                catch ME
                    disp(sprintf('error in writing, trying again | frame number: %i',i)); %#ok<DSPS>
                end
                t2 = toc;
                
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
    
end

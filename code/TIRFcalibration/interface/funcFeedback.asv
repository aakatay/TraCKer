function [fdbck] = funcFeedback(cfg,fdbck,fcall)
% feedback connections

    %% INPUTS
    % fdbck inputs : 
    % nFrst nLast runProcess syncHere syncWait toutOn isStop ssSnap ssSave 
    % fdbck inputs/outputs : 
    isSS                = fdbck.isSS;
    inSS                = fdbck.inSS;
    dispSS              = fdbck.dispSS;
    
    % worker MAT files
    MATrtWAmean         = 'signals\MATrtWAmean.mat';
    MATrtDetectThresh   = '..\signals\MATrtDetectThresh.mat';
    MATrtTraCKerPos     = '..\..\signals\MATrtTraCKerPos.mat';
    MATrtTraCKerTrace   = '..\..\signals\MATrtTraCKerTrace.mat';
    MATrtTrackSNR       = 'signals\MATrtTrackSNR.mat';
    
    %syncFrameMAT        = 'syncFrame.mat';2
    %quitToutMAT         = 'quitTout.mat';
    
    % button MAT file
    btnMAT              = 'signals\btnMAT.mat';
    
    %% associate MAT files to calling worker
    func_rtWAmean = 0;
    if strcmp(fcall,'rtWAmean')
        matFN = MATrtWAmean; 
        func_rtWAmean = 1;
    elseif strcmp(fcall,'rtDetectThresh')
        matFN = MATrtDetectThresh; 
        btnMAT = ['..\' btnMAT];
        funcIx = 2;
    elseif strcmp(fcall,'rtTraCKerPos')
        matFN = MATrtTraCKerPos;
        btnMAT = ['..\..\' btnMAT]; 
        funcIx = 3;
    elseif strcmp(fcall,'rtTraCKerTrace')
        matFN = MATrtTraCKerTrace; 
        btnMAT = ['..\..\' btnMAT];
        funcIx = 4;
    elseif strcmp(fcall,'rtTrackSNR')
        matFN = MATrtTrackSNR; 
        btnMAT = [btnMAT];
        funcIx = 5;
    end
    tloopPause = cfg.tloopPause;
    
    %% lamps (only display)
    if fdbck.runProcess
        lmpState = 1;
        save(matFN,'lmpState','-append');
    elseif fdbck.syncWait
        lmpState = -1i;
        save(matFN,'lmpState','-append');
    elseif fdbck.toutOn == 1 % timeout
        lmpState = 0;
        save(matFN,'lmpState','-append');
    elseif fdbck.syncHere
        lmpState = 1i;
        save(matFN,'lmpState','-append');
    elseif fdbck.isStop
        lmpState = -1;
        save(matFN,'lmpState','-append');
    end
    
    
    if 1 
        %% snap/save
        if dispSS % display done quit Snap/Save
            inSS = 0; dispSS = 0;
            if fdbck.ssSnap
                btnSync = 0; save(btnMAT,'btnSync','-Append'); %pushSync
            elseif fdbck.ssSave
                btnStop = 0; save(btnMAT,'btnStop','-Append'); %pushStop
            end 
            btnSnap = -1; save(btnMAT,'btnSnap','-append');  % snapping is inactive
            btnSave = -1; save(btnMAT,'btnSave','-append');  % saving is inactive
        elseif isSS
            if inSS == 0 % set button active color
                if fdbck.ssSnap, btnSnap = 1; save(btnMAT,'btnSnap','-append'); end % snapping is active
                if fdbck.ssSave, btnSave = 1; save(btnMAT,'btnSave','-append'); end % saving is active
            end 
            inSS = inSS + 1;
            if inSS > cfg.SnapNumFrames*fdbck.ssSnap + cfg.SaveNumFrames*fdbck.ssSave % stop snapping
                isSS = 0;
                dispSS = 1;
            end
        end

        %% timeout
        b_=tryLoadbtnMAT(sprintf('out=load(''%s'');',btnMAT),cfg.tTryLoop);
        if strcmp(fcall,'rtWAmean')
            if fdbck.toutOn == 1
                btnStart = -1;
                btnSync = -1;
                btnSnap = -1;
                btnSave = -1;
                btnStop = -1;    
                save(btnMAT,'btnStart', 'btnSync', 'btnSnap', 'btnSave', 'btnStop');
            elseif fdbck.toutOn == -1 % restarting
                if b_.btnStart == -1  
                    btnStart = 0;
                    save(btnMAT,'btnStart','-append');  
                end
            end
        elseif strcmp(fcall,'rtTrackSNR') 
            if fdbck.toutOn == -1 % restarted
                if b_.btnStart == 0
                    btnStart = 1;
                    save(btnMAT,'btnStart','-append');
                end
            end
        end
    end
    
    %% frame numbers
    nFrst = fdbck.nFrst;
    nLast = fdbck.nLast;
    save(matFN,'nFrst','nLast','-append');

    %% outputs:
    fdbck.isSS      = isSS;
    fdbck.inSS      = inSS;
    fdbck.dispSS    = dispSS;

% 
% 
% 
%     function funcTimeOut % paused feedback
%         while 1 % wait for btnSync
%             btnMAT_=load(btnMAT); btnSync = btnMAT_.btnSync;
%             if btnSync>0, restart = 1; noUpdate = 1; inSnap = 0; break; end
%             pause(tloopPause)
%         end
%     end
% 
% 
%     function funcFeedbackPaused % paused feedback
%         inPause = 1; inWait = 0; inWaitCounting = 0;
%         if isSnap
%             lmpState = 2; 
%         elseif isSave
%             ;
%         else
%             lmpState = 0; % timeout or sync
%         end
%         save(matFN,'lmpState','-append'); % display pause
%         while 1 % wait for btnStart
%             btnMAT_=load(btnMAT); btnStart = btnMAT_.btnStart;
%             if btnStart>0, restart = 1; noUpdate = 1; inSnap = 0; break; end
%             pause(tloopPause)
%         end
%                 
%     end % #func
    
end
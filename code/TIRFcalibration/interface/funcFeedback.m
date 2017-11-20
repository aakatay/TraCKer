function [fdbck] = funcFeedback(msgTXT,fdbck,fcall)
% feedback towards the host function
    inWait = fdbck.inWait;
    inWaitCounting = fdbck.inWaitCounting;
    inPause = fdbck.inPause;
    inSave = fdbck.inSave;
    inSaveCounting = fdbck.inSaveCounting;
    inSaveCountingIX = fdbck.inSaveCountingIX;
    inSaveCountingMAX = fdbck.inSaveCountingMAX;
    inStop = 0;

    if inWait
        if inWaitCounting
            tw = toc; % time wait
            if tw > 5 % timeout [seconds]
                funcFeedbackPaused; % check for pause and stay paused
            end
        else
            tic;
            inWaitCounting = 1;
        end
    elseif inPause
        pause(1)
        if ~exist(msgTXT), inPause = 0; end
    else
        inWait = 1;
    end
    
    fdbck.inWait = inWait;
    fdbck.inWaitCounting = inWaitCounting;
    fdbck.inPause = inPause;
    fdbck.inSave = inSave;
    fdbck.inSaveCounting = inSaveCounting;
    fdbck.inSaveCountingIX = inSaveCountingIX;
    fdbck.inStop = inStop;
    
            
            
    function funcFeedbackPaused % paused feedback
        if exist(msgTXT)
            fid = fopen(msgTXT,'r');
            ln = fgetl(fid);
            fclose(fid)
            if isempty(ln) % pause
                fTXT = sprintf(fid,'msg-paused_ %s.txt',fcall);
                fid = fopen(fTXT,'r');
                fclose(fid)
                inPause = 1; inWait = 0; inWaitCounting = 0;
                
                
            elseif strcmp(ln,'SAVE') % save
                if inSave
                    if inSaveCounting
                        inSaveCountingIX = inSaveCountingIX + 1; 
                        if inSaveCountingIX >= inSaveCountingMAX % all saved
                            fTXT = sprintf(fid,'msg-stopped_ %s.txt',fcall);
                            fid = fopen(fTXT,'r');
                            fclose(fid);        
                        end
                    else
                        inSaveCounting = 1;
                    end
                else
                    inSave = 1;
                end
                
                
            elseif strcmp(ln,'STOP') % stop
                inStop = 1;
                fTXT = sprintf(fid,'msg-stopped_ %s.txt',fcall);
                fid = fopen(fTXT,'r');
                fclose(fid)                
            end
                
            cc=3;
        end % # exist
    end % #func
    
end
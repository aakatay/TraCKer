function dbgUpdateCommunication(cfg)
% debug for rtTrackSNR.m \ updateCommunication
    % MAT files assoc. to parallel workers
    MATrtWAmean         = 'MATrtWAmean.mat';
    MATrtDetectThresh   = 'MATrtDetectThresh.mat';
    MATrtTraCKerPos     = 'MATrtTraCKerPos.mat';
    MATrtTraCKerTrace   = 'MATrtTraCKerTrace.mat';
    MATrtTraCKSNR       = 'MATrtTraCKSNR.mat';
    MATcell = {MATrtWAmean,MATrtDetectThresh,MATrtTraCKerPos,MATrtTraCKerTrace,MATrtTraCKSNR};
    
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
    
    isQuit = 0;
    
        %btnSave = 1; save(btnMAT,'btnSave','-Append');
        %btnSnap = 1; save(btnMAT,'btnSnap','-Append');
        
    while 1 % wait for startBtn
        if isequal(get(btnStart0,'BackgroundColor'), cfg.btnColPress)
            break; 
        end
        pause(0.1)
    end
    
    while 1
        if isQuit, return; end
        checkSync; % if btnSync
        checkStop; % if btnStop
        checkPostSnap; % after btnSnap
        for i = 1:5 % read inputs from five parallel workers
            matFN = MATcell{i};
            matfn = tryLoad(sprintf('out=load(''%s'');',matFN),cfg.tTryLoop); % lmpState n
            
            lmp = cfg.lmps{i};
            if matfn.lmpState==0
                set(lmp,'Color',cfg.lmpColWait); 
            elseif matfn.lmpState==2 % pause after snap 
                set(lmp,'Color',cfg.lmpColPostSnap); 
            end
            
            set(cfg.frmnoFrst{i},'Text',sprintf('%04i',matfn.nFrst)); % update text labels
            set(cfg.frmnoLast{i},'Text',sprintf('%04i',matfn.nLast));
        end
        pause(0.01)
    end


    function checkSync
        % check if btnSync
        if ~isequal(get(btnSync0,'BackgroundColor'), cfg.btnColPress), return; end
        % check if all workers paused 
        if ~isequal(get(lmp1,'Color'), cfg.lmpColWait), return; end
%         if ~isequal(get(lmp2,'Color'), cfg.lmpColWait), return; end
%         if ~isequal(get(lmp3,'Color'), cfg.lmpColWait), return; end
%         if ~isequal(get(lmp4,'Color'), cfg.lmpColWait), return; end
%         if ~isequal(get(lmp5,'Color'), cfg.lmpColWait), return; end
        % => all paused
        set(btnSync0,'BackgroundColor',cfg.btnColDefault);
        btnSync = -1; save(btnMAT,'btnSync','-Append');
        set(btnStart0,'BackgroundColor', cfg.btnColPress); % unPause
        btnStart = 1; save(btnMAT,'btnStart','-Append');
    end


    function checkPostSnap % check all workers processed SnapNumFrames 
        % check if btnStop
        if ~isequal(get(btnSnap0,'BackgroundColor'), cfg.btnColPress), return; end
        % check if all workers stopped
        if ~isequal(get(lmp1,'Color'), cfg.lmpColSnap), return; end
%         if ~isequal(get(lmp2,'Color'), cfg.lmpColStop), return; end
%         if ~isequal(get(lmp3,'Color'), cfg.lmpColStop), return; end
%         if ~isequal(get(lmp4,'Color'), cfg.lmpColStop), return; end
%         if ~isequal(get(lmp5,'Color'), cfg.lmpColStop), return; end
        % => all stopped
        set(btnSnap0,'BackgroundColor',cfg.btnColDefault);
        btnSnap = -1; save(btnMAT,'btnSnap','-Append');
        set(btnStart0,'BackgroundColor', cfg.btnColPress); % unPause
        btnStart = 1; save(btnMAT,'btnStart','-Append');
    end

    function checkStop
        % check if btnStop
        if ~isequal(get(btnStop0,'BackgroundColor'), cfg.btnColPress), return; end
        % check if all workers stopped
        if ~isequal(get(lmp1,'Color'), cfg.lmpColStop), return; end
%         if ~isequal(get(lmp2,'Color'), cfg.lmpColStop), return; end
%         if ~isequal(get(lmp3,'Color'), cfg.lmpColStop), return; end
%         if ~isequal(get(lmp4,'Color'), cfg.lmpColStop), return; end
%         if ~isequal(get(lmp5,'Color'), cfg.lmpColStop), return; end
        % => all stopped
        set(btnStop0,'BackgroundColor',cfg.btnColDefault);
        btnStop = -1; save(btnMAT,'btnStop','-Append');
        isQuit = 1;
    end
        
        
end

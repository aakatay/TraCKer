% configures parameter for realtime analysis
% remove the lines: IMG(90:94,1:4)=IMG(95:99,1:4); % ======remove
% App Designer: WEB: https://www.mathworks.com/help/matlab/components-in-app-designer.html
% image name: cy3_4nmBckgrnd_250ms_396pos_tilt8_shft3.5_[label]_0001.tif
function runRTanalysis
    F = findall(0,'type','figure'); delete(F);
    deleteALLparallelPools
    % input param
    ndigit = 4;
    bc = 4; % [nM] background concentration (for Coeff detection)
    acqTime = 0.2;% [s]
    numFrm2Save = 200; % # frames to save
    
    fclose('all');
    
    % configRT.mat
    cfg.pwd = pwd;
    cfg.msgTXT = 'msgTXT.txt';
    cfg.outDIR = 'waSeq\';
    cfg.w = [];
    cfg.h = [];
    cfg.bc = bc; % [nM] background concentration (for Coeff detection)
    cfg.waWin = 15; % walkiong average window length
    cfg.acqTime = acqTime; % [s]
    cfg.ndigit = ndigit; % # of digits for sequence number
    cfg.sptReAppearTime = 2; % [frames]
    cfg.sptJmpForTracing = 1; % [px]
    cfg.stdWin = 10; % number of frames to calc. std
    cfg.wsz = 5; % window size
    cfg.wszTracker = 5; % window size
    cfg.inSaveCountingMAX = numFrm2Save; % # frames to save
    cfg.cropTXT = [];
    cfg.cropTXT = '125X100Y50x50';
    
    
    
    %% function names
    funName = {'waSeq' 'rtDetectThresh' 'rtTraCKerPos' 'rtTraCKerTrace' 'rtTrackSNR'};
    
    fname0MAT = 'fname0.mat';
    if exist(fname0MAT), delete(fname0MAT); end
    
    %% messaging TXT file
    msgTXT = cfg.msgTXT; 
    if exist(msgTXT), delete(msgTXT); end
    

    %% buttons
    fig= uifigure('Position',[20 600 300 300]);

    x1 = 20;
    dy = 30;
    y1 = [0:3]*dy+10;
    
    w = 50;
    h1 = 20;
    btn1 = uibutton(fig,'Position',[x1 y1(1) w h1],'Text','Start','ButtonPushedFcn', @(btn1,event) runStart);
    btn2 = uibutton(fig,'Position',[x1 y1(2) w h1],'Text','Pause','ButtonPushedFcn', @(btn2,event) runPause);
    btn3 = uibutton(fig,'Position',[x1 y1(3) w h1],'Text','Save','ButtonPushedFcn', @(btn3,event) runSave);
    btn4 = uibutton(fig,'Position',[x1 y1(4) w h1],'Text','Stop','ButtonPushedFcn', @(btn4,event) runStop);

    %% text input
    x2 = w+x1+10;
    y2 = y1(1);
    w2 = 120;
    h2 = y1(end)-y1(1);
    hp = uipanel('Parent',fig,'Title','lens config','FontSize',10,...
             'BackgroundColor','white','Units','Pixels',...
             'Position',[x2 y2 w2 h1+h2]);
         
    w = 23;
	dx = w+5;
	x3 = [0:3]*dx+x2+5;
    y3 = [0:3]*dy+y2+5;
    h3 = 13;
    % rows
    labelr1 = uilabel(fig,'Position',[x3(1) y3(2) w h3],'Text','L1');
    labelr2 = uilabel(fig,'Position',[x3(1) y3(1) w h3],'Text','L2');
    % columns
    labelc1 = uilabel(fig,'Position',[x3(2) y3(3) w h3],'Text','Tilt');
    labelc2 = uilabel(fig,'Position',[x3(3) y3(3) w h3],'Text','Shft');
    labelc2 = uilabel(fig,'Position',[x3(4) y3(3) w h3],'Text','Dist');
    
    L1tilt = uieditfield(fig,'numeric','Position',[x3(2) y3(2) w h3]);
    L1shft = uieditfield(fig,'numeric','Position',[x3(3) y3(2) w h3]);
    L1dist = uieditfield(fig,'numeric','Position',[x3(4) y3(2) w h3]);
    L2tilt = uieditfield(fig,'numeric','Position',[x3(2) y3(1) w h3]);
    L2shft = uieditfield(fig,'numeric','Position',[x3(3) y3(1) w h3]);
    L2dist = uieditfield(fig,'numeric','Position',[x3(4) y3(1) w h3]);
    
    cfgLens = [L1tilt.Value L1shft.Value L1dist.Value L2tilt.Value L2shft.Value L2dist.Value];
    
    %% lamps
    x1 = 20;
    x2 = 60;
    y1 = [0:4]*30+140;

    lmp1 = uilamp(fig,'Position',[x1 y1(1) 20 20],'Color',[1 0 0]);
    lmp2 = uilamp(fig,'Position',[x1 y1(2) 20 20],'Color',[1 0 0]);
    lmp3 = uilamp(fig,'Position',[x1 y1(3) 20 20],'Color',[1 0 0]);
    lmp4 = uilamp(fig,'Position',[x1 y1(4) 20 20],'Color',[1 0 0]);
    lmp5 = uilamp(fig,'Position',[x1 y1(5) 20 20],'Color',[1 0 0]);


    label1 = uilabel(fig,'Position',[x2 y1(1) 100 15],'Text','WAmean');
    label2 = uilabel(fig,'Position',[x2 y1(2) 100 15],'Text','detectThresh');
    label3 = uilabel(fig,'Position',[x2 y1(3) 100 15],'Text','TraCKerPos');
    label4 = uilabel(fig,'Position',[x2 y1(4) 100 15],'Text','TraCKerTrace');
    label5 = uilabel(fig,'Position',[x2 y1(5) 100 15],'Text','trackSNR');

    
    
    function runStop
        fid = fopen(msgTXT,'w');
        fwrite(fid,'STOP');
        fclose(fid)
        
        %delete(msgTXT)
    end 
    
    function runSave
        save cfgLens cfgLens
        if ~exist(msgTXT) % pause
            fid = fopen(msgTXT,'w');
            fclose(fid)
            checkPause
        end
        fid = fopen(msgTXT,'w');
        fwrite(fid,'SAVE');
        fclose(fid)
    end 
    
    function checkPause
        while(1) % until all paused
            % detect pausing functions
            FN = rdir('msg-paused_*.txt');
            funPause = {};
            for i = 1:numel(fn)
                fn = FN(i).name;
                fn = fn(1:end-4); % remove '.txt'
                funPause = [funPause {sscanf(fn,'msg-paused_%s')}];
            end

            % switch lamps
            lmp = {lmp1 lmp2 lmp3 lmp4 lmp5};
            paused = 0;
            for i = 1:numel(funPause)
                fnpause = ismember(funName,funPause{i});
                if ~isempty(fnpause)
                    set(lmp{fnpause},'Color',[1 1 0]); % yellow
                end
            end

            if numel(funPause) == numel(funName) % all paused
                fn = rdir('msg-paused_*.txt'); % delete all paused messages
                for i = 1:numel(fn)
                    delete(fn(i).name);
                end
                break;
            end
            
            pause(0.1)
        end
    end
    
    
    function runPause
        
        fn = rdir('msg-paused_*.txt'); % delete all paused messages
        for i = 1:numel(fn)
            delete(fn(i).name);
        end
        
        if exist(msgTXT) % continue
            delete(msgTXT)
        else % pause
            fid = fopen(msgTXT,'w');
            fclose(fid)
        end
        checkPause;
    end 


    function runStart
        %% remove directory
        if exist(cfg.outDIR), rmdir(cfg.outDIR, 's'); end
        if exist(fname0MAT),delete(fname0MAT); end
        p = gcp();
        
        %% find crop size
        cropSize;
        % cfg.crop --> [xCr yCr szXcr szYcr]
        % reads the position of crop from the filename
        
        %% assign filename 
        while (1) % wait till first file
            fn=dir('*.tif');
            if isempty(fn), pause(0.1);continue; end
            fname0_ = fn(1).name;
            fname0 = fname0_(1:end-4); % remove .tif
            fname0 = fname0(1:end-ndigit); % remove seq number
            save(fname0MAT,'fname0');
            break;
        end
        
        nCh = numel(fname0);
        j=1;
        for i = 1: nCh
            if fname0(i) == '_'
                last_(j) = i;
                j=j+1;
            end
        end
        last_1 = last_(end-1);
        last_2 = last_(end);

        if ~exist('last_')
            label = [];
        else
            label= fname0(last_1+1:last_2-1);
        end
             
        %% save 'cfgRT.mat'
        cfg.label = label;  
        iminf = imfinfo(fname0_);
        if isempty(cfg.crop)
            cfg.w = iminf.Width;
            cfg.h = iminf.Height;
        else
            cfg.w = cfg.crop(3);
            cfg.h = cfg.crop(4);
        end
        
        save cfgRT cfg;    
        
        %% verify running codes
        p = gcp();f1=parfeval(p,@rtWAmean,0);
        while (1), if exist('waSeq'), break; end;end
        set(lmp1,'Color',[0 1 0]);
        p = gcp();f2=parfeval(p,@rtDetectThresh,0);
        while (1), if exist('waSeq\tracker'), break; end; end
        set(lmp2,'Color',[0 1 0]);
        p = gcp();f3=parfeval(p,@rtTraCKerPos,0);
        while (1), if ~isempty(rdir('waSeq\tracker\posData-coeff*.mat')), break; end; end
        set(lmp3,'Color',[0 1 0]);
        p = gcp();f4=parfeval(p,@rtTraCKerTrace,0);
        while (1), if ~isempty(rdir('waSeq\tracker\cfgRT')), break; end; end
        set(lmp4,'Color',[0 1 0]);
        p = gcp();f5=parfeval(p,@rtTrackSNR,0);
        while (1), if ~isempty(rdir('waSeq\tracker\cfgRT\SMdata.mat')), break; end; end
        set(lmp5,'Color',[0 1 0]);

        if 0
            deleteALLparallelPools
        end
    end

    function cropSize
        % cfg.cropTXT --> cfg.Crop ([xCr yCr szXcr szYcr])
        % reads the position of crop from the filename
        cropTXT = cfg.cropTXT;
        if isempty(cropTXT), return; end
        
        nCh = numel(cropTXT);
        
        nmCrop = cropTXT;
        ixX = find(nmCrop=='X');
        ixY = find(nmCrop=='Y');
        ixx = find(nmCrop=='x');

        xCr = str2num(nmCrop(1:ixX-1));
        yCr = str2num(nmCrop(ixX+1:ixY-1));
        szXcr = str2num(nmCrop(ixY+1:ixx-1));
        szYcr = str2num(nmCrop(ixx+1:end));
        cfg.crop = [xCr yCr szXcr szYcr];
    end
end

                
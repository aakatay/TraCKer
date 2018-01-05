% configures parameter for realtime analysis
% remove the lines: IMG(90:94,1:4)=IMG(95:99,1:4); % ======remove
% App Designer: WEB: https://www.mathworks.com/help/matlab/components-in-app-designer.html
% image name: cy3_4nmBckgrnd_250ms_396pos_tilt8_shft3.5_[label]_0001.tif
function runRTanalysis
    isdbgAcquisition = 1;
    F = findall(0,'type','figure'); delete(F);
    %deleteALLparallelPools
    fclose('all');
    gcp
    % input param
    ndigit = 4;
    bc = 4; % [nM] background concentration (for Coeff detection)
    acqTime = 1;% [s]
    acqTime = .4;% [s]
    numFrm2Save = 200; % # frames to save
    numFrm2Snap = 20;
    
    
    % configRT.mat
    
    cfg.camPixelSize = 16e3; % [nm] camera pixel size
    cfg.imagingMag = 100*1.5; % magnification (Kural lab)
    cfg.imagingMag = 60; % magnification (Fishel lab)
    cfg.imagingPixelSize = cfg.camPixelSize/cfg.imagingMag;
    cfg.pwd = pwd;
    cfg.msgTXT = 'msgTXT.txt';
    cfg.outDIR = 'waSeq\';
    cfg.acqDIR = 'acq\';
    cfg.snapSaveDIR = 'snapSaveOUT\';
    cfg.signalDIR = 'signals\';
    cfg.w = [];
    cfg.h = [];
    cfg.bc = bc; % [nM] background concentration (for Coeff detection)
    cfg.waWin = 15; % walkiong average window length
    cfg.acqTime = acqTime; % [s]
    cfg.ndigit = ndigit; % # of digits for sequence number
    cfg.sptReAppearTime = 2; % [frames]
    cfg.sptJmpForTracing = 1; % [px]
    cfg.stdWin = 10; % number of frames to calc. std
    cfg.wsz = 6; % [even number] window size
    cfg.wsz2 = 2; % [even number] intensity window size
    cfg.wszTracker = 5; % window size
    cfg.numFrm2Save = numFrm2Save; % # frames to save
    cfg.numFrm2Snap = numFrm2Snap; % # frames to snap
    cfg.cropTXT = [];
    %cfg.cropTXT = '125X100Y50x50';
    %cfg.cropTXT = '0X0Y512x512';
    cfg.cropTXT = '0X0Y260x260';
    cfg.dispMag = 4;
    cfg.logWA = [pwd '\logData\logWA.txt'];
    cfg.logThresh = [pwd '\logData\logThresh.txt'];
    cfg.logPos = [pwd '\logData\logPos.txt'];
    cfg.logTrace = [pwd '\logData\logTrace.txt'];
    cfg.logSNR = [pwd '\logData\logSNR.txt'];
    cfg.SNRcolorThresh = [1.5 7 12];
    cfg.mxNumLocalization = 2*1e3; % need to change ->change the code rtTraCKerTrace.m
    cfg.scrnSzIn = [1600 900];
    cfg.isTlog = 1; % keeps times in TXT files in logData\
    cfg.tloopPause = 0.01; % [sec]
    cfg.timeOut = 10; % [sec]
    cfg.timeOut = inf; % [sec]
    cfg.tTryLoop = 0.01;  % [sec]
    cfg.SnapNumFrames = 30;  % [frames]
    cfg.SaveNumFrames = 60;  % [frames]
    cfg.isdbgAcquisition = isdbgAcquisition;
    
    % colors (buttons)
    cfg.btnColDefault   = [0.9600 0.9600 0.9600 ]; % gray
    cfg.btnColPress     = [1.0 1.0 0.0]; % yellow
    cfg.btnColActive    = [0.1 0.9 0.1]; % green
    cfg.btnColfrstData  = [0 1 1]; % light blue
    
    % colors (lamps)
    cfg.lmpColStop      = [1 0 0 ]; % red
    cfg.lmpColSyncWait  = [1 0.5 0 ]; % orange
    cfg.lmpColTimeout   = [0 0 0]; % black
    cfg.lmpColSyncHere  = [1 0 1]; % magenta
    cfg.lmpColActive    = [0 1 0]; % green
    
    %% folder structure
    if 1 
        if exist('acq'), rmdir('acq', 's'); end; mkdir('acq')
        if exist('snapSaveOUT'), rmdir('snapSaveOUT', 's'); end; mkdir('snapSaveOUT')
        if exist('signals'), rmdir('signals', 's'); end; mkdir('signals')
        if exist(cfg.outDIR), rmdir(cfg.outDIR, 's'); end
        mkdir waSeq
        mkdir waSeq\tracker
        mkdir waSeq\tracker\rtData
        while ~exist('waSeq\tracker\rtData')
            pause(0.1)
        end
    end

    %% communication
    MATrtWAmean         = 'signals\MATrtWAmean.mat';
    MATrtDetectThresh   = 'signals\MATrtDetectThresh.mat';
    MATrtTraCKerPos     = 'signals\MATrtTraCKerPos.mat';
    MATrtTraCKerTrace   = 'signals\MATrtTraCKerTrace.mat';
    MATrtTrackSNR       = 'signals\MATrtTrackSNR.mat';
    btnMAT              = 'signals\btnMAT.mat';
    
    lmpState = -1; % (-1: stop 0:pause 1:active)
    nFrst = 0;
    nLast = 0;
    btnStart = -1;
    btnSync = -1;
    btnSnap = -1;
    btnSave = -1;
    btnStop = -1;    
    
    save(btnMAT,'btnStart', 'btnSync', 'btnSnap', 'btnSave', 'btnStop');
    save(MATrtWAmean,       'lmpState', 'nFrst', 'nLast');
    save(MATrtDetectThresh, 'lmpState', 'nFrst', 'nLast');
    save(MATrtTraCKerPos,   'lmpState', 'nFrst', 'nLast');
    save(MATrtTraCKerTrace, 'lmpState', 'nFrst', 'nLast');
    save(MATrtTrackSNR,     'lmpState', 'nFrst', 'nLast');
    
    % parallel workers
    f1=[];f2=[];f3=[];f4=[];
    
    %% log
    fclose('all');
    if exist('logData')~=7
        mkdir('logData'); 
    else
        rmdir('logData', 's');
        mkdir('logData'); 
    end
    
    
    %% voronoi figure
    nlast = 0;
    figSNRvoronIMG = [];
    CM = gray(256);
    szXY = [];
    mag = [];
    if ~exist('output'), mkdir('output'); end
    SNRmovieVoronoiFN = 'output\SNRmovieVoronoi.tif'; % output
    
    %% function names
    funName = {'waSeq' 'rtDetectThresh' 'rtTraCKerPos' 'rtTraCKerTrace' 'rtTrackSNR'};
    
    fname0MAT = 'fname0.mat';
    if exist(fname0MAT), delete(fname0MAT); end
    
    %% messaging TXT file
    msgTXT = cfg.msgTXT; 
    if exist(msgTXT), delete(msgTXT); end
    

    %% buttons
    fig= uifigure(1,'Position',[20 600 300 400]);
    

    x1 = 20;
    dy = 30;
    y1 = [0:4]*dy+10;
    
    w = 50;
    h1 = 20;
    btnStart0 = uibutton(fig,'Position',[x1 y1(1) w h1],'Text','Start','ButtonPushedFcn', @(btnStart0,event) runStart);
    btnSync0 = uibutton(fig,'Position',[x1 y1(2) w h1],'Text','Sync','ButtonPushedFcn', @(btnSync0,event) runSync);
    btnSnap0 = uibutton(fig,'Position',[x1 y1(3) w h1],'Text','Snap','ButtonPushedFcn', @(btnSnap0,event) runSnap);
    btnSave0 = uibutton(fig,'Position',[x1 y1(4) w h1],'Text','Save','ButtonPushedFcn', @(btnSave0,event) runSave);
    btnStop0 = uibutton(fig,'Position',[x1 y1(5) w h1],'Text','Stop','ButtonPushedFcn', @(btnStop0,event) runStop);

    cfg.btns.btnStart0 = btnStart0;
    cfg.btns.btnSync0 = btnSync0;
    cfg.btns.btnSnap0 = btnSnap0;
    cfg.btns.btnSave0 = btnSave0;
    cfg.btns.btnStop0 = btnStop0;
    
    %% text input
    x2 = w+x1+10;
    y2 = y1(1);
    w2 = 120;
    h2 = y1(end)-y1(1);
    hT1 = h1+h2-40;
    hp = uipanel('Parent',fig,'Title','lens config','FontSize',10,...
             'BackgroundColor','white','Units','Pixels',...
             'Position',[x2 y2 w2 hT1]);
         
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
    
    cfg.lensBox = [L1tilt L1shft L1dist L2tilt L2shft L2dist];
    cfgLens = [L1tilt.Value L1shft.Value L1dist.Value L2tilt.Value L2shft.Value L2dist.Value];
    
    %% TIRF angle
    
    dh = 90;
    dxx = 10;
    hp2 = uipanel('Parent',fig,'FontSize',10,...
             'BackgroundColor','white','Units','Pixels',...
             'Position',[x2 y2+hT1 w2 40]);
         
    labelt1 = uilabel(fig,'Position',[x3(1) y3(2)+dh w*2 h3],'Text','TIRFx');
    labelt2 = uilabel(fig,'Position',[x3(1) y3(1)+dh+12 w*2 h3],'Text','ACQx');
    TIRFx0 = uieditfield(fig,'numeric','Position',[x3(3)-10 y3(2)+dh w*2.5 h3]);
    TIRFx = uieditfield(fig,'numeric','Position',[x3(3)-10 y3(1)+dh+12 w*2.5 h3]);

    %% lamps
    dxfn = 50; % frame num width
    x1 = 10;
    x2 = x1+35;
    x3 = x2+30;
    x4 = x3+40;
    y1 = [0:4]*30+170;

    % frame number (first)
    frmnoFrst1 = uilabel(fig,'Position',[x1 y1(1) 100 15],'Text','0000');
    frmnoFrst2 = uilabel(fig,'Position',[x1 y1(2) 100 15],'Text','0000');
    frmnoFrst3 = uilabel(fig,'Position',[x1 y1(3) 100 15],'Text','0000');
    frmnoFrst4 = uilabel(fig,'Position',[x1 y1(4) 100 15],'Text','0000');
    frmnoFrst5 = uilabel(fig,'Position',[x1 y1(5) 100 15],'Text','0000');
    %lamps
    lmp1 = uilamp(fig,'Position',[x2 y1(1) 20 20],'Color',[1 0 0]); % wamean
    lmp2 = uilamp(fig,'Position',[x2 y1(2) 20 20],'Color',[1 0 0]);
    lmp3 = uilamp(fig,'Position',[x2 y1(3) 20 20],'Color',[1 0 0]);
    lmp4 = uilamp(fig,'Position',[x2 y1(4) 20 20],'Color',[1 0 0]);
    lmp5 = uilamp(fig,'Position',[x2 y1(5) 20 20],'Color',[1 0 0]);
    
    % frame number (last)
    frmnoLast1 = uilabel(fig,'Position',[x3 y1(1) 100 15],'Text','0000');
    frmnoLast2 = uilabel(fig,'Position',[x3 y1(2) 100 15],'Text','0000');
    frmnoLast3 = uilabel(fig,'Position',[x3 y1(3) 100 15],'Text','0000');
    frmnoLast4 = uilabel(fig,'Position',[x3 y1(4) 100 15],'Text','0000');
    frmnoLast5 = uilabel(fig,'Position',[x3 y1(5) 100 15],'Text','0000');
    % error buttons
    dy=3;
    label1 = uibutton(fig,'Position',[x4 y1(1)-dy 100 15+2*dy],'Text','WAmean','ButtonPushedFcn', @(btnStart,event) queryWAmean);
    label2 = uibutton(fig,'Position',[x4 y1(2)-dy 100 15+2*dy],'Text','detectThresh','ButtonPushedFcn', @(btnStart,event) querydetectThresh);
    label3 = uibutton(fig,'Position',[x4 y1(3)-dy 100 15+2*dy],'Text','TraCKerPos','ButtonPushedFcn', @(btnStart,event) queryTraCKerPos);
    label4 = uibutton(fig,'Position',[x4 y1(4)-dy 100 15+2*dy],'Text','TraCKerTrace','ButtonPushedFcn', @(btnStart,event) queryTraCKerTrace);
    label5 = uilabel(fig,'Position',[x4 y1(5) 100 15],'Text','trackSNR');
    
    cfg.lmps = {lmp1, lmp2, lmp3, lmp4, lmp5};
    cfg.frmnoFrst = {frmnoFrst1, frmnoFrst2, frmnoFrst3, frmnoFrst4, frmnoFrst5};
    cfg.frmnoLast = {frmnoLast1, frmnoLast2, frmnoLast3, frmnoLast4, frmnoLast5};
    
	%dbgUpdateCommunication(cfg);
    %while 1, pause(0.1); end
    runStart
    return
    
    
%% functions ====================================================    
    % get worker errors
    function queryWAmean
        errTxt = 'no error';
        if ~isempty(f1.Error), errTxt = getReport(f1.Error); end
        msgbox(sprintf('queryWAmean\n : STATE:%s\n ERROR:%s\n ',f1.State,errTxt));
    end
    function querydetectThresh
        errTxt = 'no error'; if ~isempty(f2.Error), errTxt = getReport(f2.Error); end
        msgbox(sprintf('querydetectThresh\n : STATE:%s\n ERROR:%s\n ',f2.State,errTxt));
    end
    function queryTraCKerPos
        errTxt = 'no error'; if ~isempty(f3.Error), errTxt = getReport(f3.Error); end
        msgbox(sprintf('queryTraCKerPos\n : STATE:%s\n ERROR:%s\n ',f3.State,errTxt));
    end
    function queryTraCKerTrace
        errTxt = 'no error'; if ~isempty(f4.Error), errTxt = getReport(f4.Error); end
        msgbox(sprintf('queryTraCKerTrace\n : STATE:%s\n ERROR:%s\n ',f4.State,errTxt));
    end
    
    %buttons
    function runStop
        btnStop = 0; save(btnMAT,'btnStop','-Append');
    end 
    
    function runSave
        btnSave = 0; save(btnMAT,'btnSave','-Append');
        runSync
    end 

    function runSnap
        btnSnap = 0; save(btnMAT,'btnSnap','-Append');
        runSync
    end
    
    function runSync
        btnSync = 0; save(btnMAT,'btnSync','-Append');
    end 

    % start
    function runStart
        if strcmp(get(btnStart0,'BackgroundColor'), cfg.btnColPress), return; end % already running
        if strcmp(get(btnStart0,'BackgroundColor'), cfg.btnColActive), return; end % already running
        btnStart = 0; save(btnMAT,'btnStart','-Append');
        
        set(btnStart0,'BackgroundColor', cfg.btnColfrstData);
        runSync
        %% remove directory
        if exist(fname0MAT),delete(fname0MAT); end
        
        %% find crop size
        cropSize;
        % cfg.crop --> [xCr yCr szXcr szYcr]
        % reads the position of crop from the filename
        
        if isdbgAcquisition
            inputFN_ = dir(['input\*.tif']);
            inputFN0 = inputFN_(1).name;
            inputFN = ['input\' inputFN0];
            IMGin = imread(inputFN,1);
            inputFN2  = sprintf('acq\\%s_%04i.tif',inputFN0(1:end-4),1);
            imwrite(IMGin,inputFN2);
        end
        
        %% assign filename 
        while (1) % wait till first file
            fn=dir('acq\*.tif');
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
        iminf = imfinfo(['acq\' fname0_]);
        if isempty(cfg.crop)
            cfg.w = iminf.Width;
            cfg.h = iminf.Height;
        else
            cfg.w = cfg.crop(3);
            cfg.h = cfg.crop(4);
        end
        
        cfgSave = cfg;
        cfgSave.btns = [];
        cfgSave.lmps = [];
        cfgSave.lensBox = [];
        save('cfgRT','cfgSave');    
        
        %% prepare
        %prepFigSNRdata;
        
        %% verify running codes
        p = gcp('nocreate');f1=parfeval(p,@rtWAmean,0,cfgSave);
        p = gcp('nocreate');f2=parfeval(p,@rtDetectThresh,0,cfgSave);
        p = gcp('nocreate');f3=parfeval(p,@rtTraCKerPos,0,cfgSave);
        p = gcp('nocreate');f4=parfeval(p,@rtTraCKerTrace,0,cfgSave);
        %return; 
        rtTrackSNR(cfg)
        %p = gcp('nocreate');f5=parfevalOnAll(p,@rtTrackSNR,0);
        ccc=3;
%         return
%         while (1), if ~isempty(rdir('waSeq\*.tif')), break; end; pause(0.5);end
%         set(lmp1,'Color',[0 1 0]);
%         while (1), if ~isempty(rdir('waSeq\tracker\Coeff*.mat')), break; end; pause(0.5); end
%         set(lmp2,'Color',[0 1 0]);
%         while (1), if ~isempty(rdir('waSeq\tracker\posData-coeff*.mat')), break; end; pause(0.5);end
%         set(lmp3,'Color',[0 1 0]);
%         while (1), if ~isempty(rdir('waSeq\tracker\rtData\traceData_*.mat')), break; end; pause(0.5); end
%         set(lmp4,'Color',[0 1 0]);
% 
%         if 0
%             deleteALLparallelPools
%         end
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


















% junk vvvvvvvvvvvvvvvvvvvvvvvvvv

    function prepFigSNRdata
        return;
        w = cfg.w;
        h = cfg.h;
        szXY = [w h];
        mag = cfg.dispMag; % display size
        
        %% SNR image figure
        [mag, pos, szx, szy ] = calcMaxMag(zeros(szXY),mag,cfg.scrnSzIn);
        szXYmag = [szx, szy];
        szYXmag = fliplr(szXYmag);

        tit = 'SNR voronoi image';
        pos(1) = pos(1) - 800;
        pos(1) = pos(1)- 700;    
        figSNRvoronIMG = figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos/2 szXYmag(1) szXYmag(2)]);
        axeSNRvoron = axes('Parent',figSNRvoronIMG,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren','XLim',0.5+[0 szXYmag(1)],'YLim',0.5+[0 szXYmag(2)]);

    end

    function dispSNRvoronoiIMG
        SNRdataFN = [pwd '\waSeq\tracker\rtData\SNRdata.mat'];
        SNRdata = load(SNRdataFN);
        ixFrm = SNRdata.ixFrm;
        XYS = SNRdata.XYS;
        n = nlast;
        xys = XYS(ixFrm(n-1)+1:ixFrm(n),:);
        if size(xys,1)<5 % write a blank image
            SNRmovVoronoi = zeros(szXY*mag);
            %imwrite(SNRmovVoronoi,SNRmovieVoronoiFN,'WriteMode','append','Compression', 'none') 
        else
            snrDisp = xys(:,3);

            % SNR intensity range scales
            if max(snrDisp)<1.5 
                sscale = 1.5;
                EdgeColorSel = 1; % red
            elseif max(snrDisp)<5
                sscale = 5;
                EdgeColorSel = 2; % green 
            else % max(snrDisp)<10
                sscale = 10;
                EdgeColorSel = 3; % blue
            end
            [SNRmovVoronoi,~] = getVoronoinImg(figSNRvoronIMG,xys(:,1:2),snrDisp/sscale,szXY,mag,CM,EdgeColorSel);
            imwrite(SNRmovVoronoi,SNRmovieVoronoiFN,'WriteMode','append','Compression', 'none') 
        end
    end

    function bool = waitSNRdata 
        bool = 0;
        SNRdataFN = [pwd '\waSeq\tracker\rtData\SNRdata.mat'];
        SNRdata_ = dir(SNRdataFN);
        if isempty(SNRdata_), return; end
        SNRdata=load(SNRdataFN, 'ixFrm');
        ixFrm = SNRdata.ixFrm;
        nInput = numel(ixFrm);
        if nlast == nInput
            return;
        else
            nlast = nInput;
            snrFrames = [snrFrames nlast];
            save(snrFrm,'snrFrames');
            bool = 1; % continue
        end    
        
    end
end

                
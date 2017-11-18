clear; close all;
F = findall(0,'type','figure','tag','TMWWaitbar'); delete(F)
% searches subfolders for trace data and display combined data
% run in cell folder with subfolders such as:
% tirf001_crop4\32X25Y44x40\_000-coeff762
% make sure there is only one data set per folder

intFilter = 0; % intensity filter
isTrace=2; % ============================ select overlay data =====================================
%1:spots 2:traces 3:trace average 4:trace recruit average (localized traces) 

%% load files
% image folder
isWrtLongTraces = 0;

imgsFolder = brdir('imgs'); % recursive folder search
fn = dir('img*.tif');
imgFile = [];
if isempty(fn)
    fn = dir([imgsFolder '\img*.tif']);
    if isempty(fn)
        %warning('no image file exists')
    else
        imgFile = [imgsFolder '\' fn(1).name];
    end
else
    imgFile = fn(1).name;
end

% parameters

isTraceCombination = 0;
resMag = 1;
genImg = 0;
genOverlayImg = 0;
isDispRecruitment = 1;
%if genImg, genOverlayImg = 0; end;
isCountInc = 1; % every recruitment counted as one
highResImg = 'highResImg.tif';
% pixelation
frmBinDetect = 1; % param used in tracker
pxBin = 0.25; % pixel size multiplier

fnXYZtraceData0 = rdir('**\*traceData0-coeff*');
fnXYZ = rdir('**\*xyzDataGaus-coeff*');
% combine trace files
%fnXYZtraceData0 = flipud(fnXYZtraceData0);
%fnXYZ = flipud(fnXYZ);
if numel(fnXYZtraceData0)==1
    isCombineDataSets = 0;
end

TraceX_ = [];
TraceY_ = [];
TraceINT_ = [];
frmNoTrace_ = [];
trInf_ = [];
frames = 0;
tracePos = 0;

for i=1:numel(fnXYZtraceData0) % combine data sets
    load(fnXYZtraceData0(i).name);
    load(fnXYZ(i).name,'ixSptFrm');
    % frames cropping
%cfg.img.lastFrm = 30000;
    frmCrp(i,:) = [cfg.img.frstFrm cfg.img.lastFrm];
    disp(sprintf('crop#%i: frames %i-%i',i,frmCrp(i,1),frmCrp(i,2) ));
    TraceX_ = [TraceX_ TraceX];
    TraceY_ = [TraceY_ TraceY];
    TraceINT_ = [TraceINT_ TraceINT];
    frmNoTrace_ = [frmNoTrace_ frames+frmNoTrace];
    trInf(:,1) = trInf(:,1) + frames;
    trInf(:,3) = trInf(:,3) + tracePos;
    trInf_ = vertcat(trInf_,trInf);
    frames = frames+frmCrp(i,2)-frmCrp(i,1)+1;
    tracePos = tracePos+numel(TraceX);
    Frames(i) = numel(ixSptFrm)-1;
    disp(sprintf('data folder #: %i, #frames:%i, #traces:%i',i,Frames(i),size(TraceX,2)))

    % folder names
    slashPos = strfind(fnXYZtraceData0(i).name,'\');
    folderName = fnXYZtraceData0(i).name;
    folderNames(i) = {folderName(1:slashPos(end))};
end
TraceX = TraceX_;
TraceY = TraceY_;
TraceINT = TraceINT_;
frmNoTrace = frmNoTrace_;
trInf = trInf_;

%load('linearData'); % CMEanalysis
%load('cfg')
%frmCrp(1,:) = [cfg.img.frstFrm cfg.img.lastFrm];
%frames = frames+frmCrp(1,2)-frmCrp(1,1)+1;
% CD = cd;
% CoeffNormAll = [];
% for i=1:numel(fnXYZtraceData0) % join coeff values
%     cd(cell2mat(folderNames(i)));
% %    load CoeffFit;
% %    Coeff = CoeffFit(1);
%     %load CoeffNormMulSmth;
% %    CoeffNorm = CoeffFit;
% %    CoeffNormAll = [CoeffNormAll CoeffNorm];
%     cd(CD);
% end
   
% stats
if ~exist('.\stats'), mkdir('stats'), end;
%save('stats\CoeffNormAll','CoeffNormAll')

%plot(CoeffNormAll); title('Coefficient change'); ylabel('Coefficient'); xlabel('frames');
%mx= max(CoeffNormAll);

if 0 
    FramesLast = 0;
    for i=1:numel(fnXYZtraceData0) % join coeff values
        if i > 1, FramesLast = sum(frmCrp(1:i-1,2)); end
        hl=line([FramesLast+frmCrp(i,1) FramesLast+frmCrp(i,1)],[0 mx]);
        line([FramesLast+frmCrp(i,2) FramesLast+frmCrp(i,2)],[0 mx]);
        set(hl,'LineStyle','-.','Color',[1 0 0])
    end
end

imgFig = getframe(gcf);
imgOut = imgFig.cdata;
imwrite(imgOut,'CoeffNormAll.tiff','Compression', 'none')
close(gcf)


if 0
%% temporarily add info from files
    % read cell info
    cfg.cell = struct;
    fcell = 'info-cell.txt';
    if isempty(fcell), error('no cell info is found'); end;
    hfcell = fopen(fcell);
    cfg.cell.cellType = fscanf(hfcell,'cellType:%s\n');
    cfg.cell.platingTime = fscanf(hfcell,'platingTime:%s\n');
    cfg.cell.surfTreat = fscanf(hfcell,'surfTreat:%s\n');
    cfg.cell.chemTreat = fscanf(hfcell,'chemTreat:%s\n');
    cfg.cell.other = fscanf(hfcell,'other:%s');
    fclose(hfcell);
    % read exp info
    cfg.exp = struct;
    fexp = 'info-exp.txt';
    if isempty(fexp), error('no exp info is found'); end;
    hfexp = fopen(fexp);
    cfg.exp.date = fscanf(hfexp,'date:%s\n');
    cfg.exp.time = fscanf(hfexp,'time:%s\n');
    cfg.exp.other = fscanf(hfexp,'other:%s');
    fclose(hfexp);
    % read bleach info
    cfg.bleach = struct;
    fbleach = 'info-bleach.txt';
    if isempty(fbleach), error('no bleach info is found'); end;
    hfbleach = fopen(fbleach);
    cfg.bleach.imgName = fscanf(hfbleach,'pre-image name:%s\n');
    cfg.bleach.acqTime = fscanf(hfbleach,'acq. time:%s\n');
    i = 1;
    while ~feof(hfbleach)
        cfg.bleach.bleachData{i,1} = fscanf(hfbleach,'ND filter position:%s\n',1);
        cfg.bleach.bleachData{i,2} = fscanf(hfbleach,'bleach time:%s\n',1);
        cfg.bleach.bleachData{i,3} = fscanf(hfbleach,'TIRF position:%s\n',1);
        i = i + 1;
    end
    fclose(hfbleach);
    % read acquisition info
    cfg.acq = struct;
    facq = 'info-acq.txt';
    if isempty(facq), error('no acquisition info is found'); end;
    hfacq = fopen(facq);
    cfg.acq.ROIcam = fscanf(hfacq,'camROI: %s\n');
    cfg.acq.ROInd = fscanf(hfacq,'ndROI:%s\n');
    cfg.acq.TIRFstagePosEPI = fscanf(hfacq,'EPIstagePos:%s\n');
    cfg.acq.TIRFstagePosTIRF = fscanf(hfacq,'TIRFstagePos:%s\n');
    cfg.acq.tempBeforePositioning = fscanf(hfacq,'1-tempBeforePositioning:%s\n');
    cfg.acq.tempChangeInTenMin = fscanf(hfacq,'2-tempChangeInTenMin:%s\n');
    cfg.acq.power = fscanf(hfacq,'power:%s\n');
    cfg.acq.gain = fscanf(hfacq,'gain:%s\n');
    cfg.acq.NDfilt = fscanf(hfacq,'NDfilt:%s\n');
    cfg.acq.pixelCAM = fscanf(hfacq,'pixelCAM:%s\n');
    cfg.acq.magObj = fscanf(hfacq,'magObj:%s\n');
    cfg.acq.magEx	 = fscanf(hfacq,'ExMag:%s\n');
    cfg.acq.other = fscanf(hfacq,'other:%s');
    fclose(hfacq);
end
pixelCAM                = str2num(cfg.acq.pixelCAM);
magObj                  = str2num(cfg.acq.magObj);
magEx                   = str2num(cfg.acq.magEx);

pxSz =  round(pixelCAM/magObj/magEx*10000)/10; %nm
%% load filename
%slashFind=strfind(fnXYZtraceData0(end).name,'\');
%folder = fnXYZtraceData0(end).name;
% if isempty(slashFind)
%     load('fname.mat')
% else
%     load([folder(1:slashFind(end)-1) '\fname.mat'])
% end


if ~exist('cfg'), cfg = struct; end;
cfg.overlay = struct;
if ~isfield(cfg,'report'), cfg.report = struct; end;

%frames = 13107;
if 0
    % exposure time
    msPos2 = strfind(fname,'ms');
    fname_ = fname(1:msPos2);
    msPos1 = strfind(fname_,'_');
    msPos1 = msPos1(end)+1;
    acqTime = str2num(fname(msPos1:msPos2-1))/1000; % sec
end
acqTime = cfg.acq.acqTime;
maxNumFrames_time = 120; % 2min
if 120/acqTime > frames
    maxNumFrames_time = frames*acqTime;
end

% update minTrLen
binFrame = cfg.img.binFrame;
minTrLen = 2;
maxTrLenTime = 200; % ms
maxTrLen = floor(maxTrLenTime/acqTime);

cfg.overlay.isCountInc = isCountInc;
%% data cropping 12X52Y47x47
isDataCropped=0; % = 1 if needs to be cropped
isNDcrop = 0; % ow imageJ crop
cfg.overlay.crop = struct;
cfg.overlay.crop.isNDcrop = isNDcrop;


a=rdir('**\*fname.mat');
load(a(1).name)


nCh = numel(fname);
last_= strfind(fname,'_');
%last = last_(end-2:end-1);
%last = last_(end-1:end);
if numel(last_) > 1
    %last = last_(end-2:end-1);
    %last = last_(end-1:end);
    nmCrop = fname(last_(end-2)+1:last_(end-1)-1);
    %nmCrop = fname(last_(end-1)+1:last_(end)-1);
    ixX = find(nmCrop=='X');
    ixY = find(nmCrop=='Y');
    ixx = find(nmCrop=='x');
end

% ixX = find(nmCrop=='x');
% ixY = find(nmCrop=='y');
% ixx = find(nmCrop=='X');

if ~isempty(ixX) && ~isempty(ixY) && ~isempty(ixx) % image needs to  be cropped
    isDataCropped = 1;
    xCr = str2num(nmCrop(1:ixX-1));
    yCr = str2num(nmCrop(ixX+1:ixY-1));
    szXcr = str2num(nmCrop(ixY+1:ixx-1));
    szYcr = str2num(nmCrop(ixx+1:end));
    cfg.overlay.crop.xCr_X_yCr_Y_szXcr_x_szYcr = [xCr yCr szXcr szYcr]; 
end
cfg.overlay.crop.isDataCropped = isDataCropped;

%     xParse = sscanf(cfg.acq.ROIcam,'%dX%dY%dx%d');
%     dx = xParse(3);
%     dy = xParse(4);
%     a = zeros(dx,dy);
if isempty(imgFile)
    %a = zeros(125,125);
    a = zeros(szYcr,szXcr);
else
    a=imread(imgFile);
    if isDataCropped, a= a(yCr+1:yCr+szYcr,xCr+1:xCr+szXcr); end;
end

%a = zeros(101,55);
mag = 15;
[mag, pos, m, n ] = calcMaxMag(a,mag);

memSize = size(a,2)*size(a,1)*frames*64;
isOverSize = 0;
if memSize/1e6 > 3000
    isOverSize = 1;
end
q = 0;
q_ = 0;

traceORspot{1} = 'Spot';
traceORspot{2} = 'Trace';
traceORspot{3} = 'TraceAv'; 
traceORspot{4} = 'TraceRecAv';

    
%% trace - spot loop
if genOverlayImg
    tit = 'reconstructed image';
    CM = genColorMap('jet',frames); 
    figImg=figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos/2 m n]);
    axis image
    axe=axes('Parent',figImg,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off');

    hImg = imagesc(a);
    axis image
    scX2 = size(a,2)*19/20;
    scY = size(a,1)*19/20;
    scLen = 200/pxSz; % scale for 200nm
    scX = scX2-scLen;
    scY2 = scY;

    hScale = line([scX2 scX],[scY scY2],'LineWidth',2,'Color',[1 1 1]);
    hTextScale= text(scX,scY-2,'200nm','Color',[1 1 1]);
end

if isTrace>1
    traceX = sparse(size(TraceX,1),size(TraceX,2));
    traceY = sparse(size(TraceY,1),size(TraceY,2));

    
    %trace length
    figure;
    [bins, binPos] = hist(trInf(:,2),max(trInf(:,2)));
    bins = bins  ./ sum(bins); bar(binPos,bins);
    set(gca,'Color',[1 1 1]); title(sprintf('trace length - minTraceLen:%i',minTrLen)); xlabel('# of frames'); ylabel('# of traces (normalized)');
    hl=line([minTrLen,minTrLen],[0,max(bins)]); set(hl,'Color',[1 0 0 ]);
    imgFig = getframe(gcf); imgOut = imgFig.cdata;
    imwrite(imgOut,'stats\trace_lengths.tif');
    %maximize
    
    %pause;
    close(gcf)
    pause(0.5)
    
    if minTrLen > 2 
        trInf = trInf(trInf(:,2)>=minTrLen,:); % select long traces
    end
    trInf2 = trInf;
    cfg.overlay.minTrLen = minTrLen;
    

%            traceRecX = sparse(size(TraceX,1),size(TraceX,2));
%            traceRecY = zeros(size(TraceY,1),size(TraceY,2));
    nFail = 0;
   
    %disp(sprintf('trace dispersion ellimination: %.02f percent is discarded',(size(trInf,1)-size(trInf2,1))/size(trInf,1)*100));
    cfg.report.disperseTraceElimRate = (size(trInf,1)-size(trInf2,1))/size(trInf,1);
end



switch isTrace
    case 1 % spots
        load(fnXYZ.name)
        frmNo = frmNoSpot;
        textDisp = 'incidences';
    case 2 % traces
        frmNo = frmNoTrace;
        textDisp = 'incidences in the traces';
    case 3 % trace average
        frmNo = frmNoTrace;
        textDisp = 'trace average';
    case 4 % trace recruit average (localized traces)
        frmNo = frmNoTrace;
        textDisp = 'trace recruit average (only localized traces)';
end    

% correct TraceInt values
TraceInt = TraceINT; % ---
TraceInt(TraceInt<0)=0; TraceInt(isnan(TraceInt))=0;


%clear X Y INT TraceX TraceY TraceINT traceX traceX traceInt traceRecX traceRecY traceRecInt;
%xx = xx-0.5; yy = yy-0.5;
fr1 = 1;
fr2 = frames;
numFrames = fr2 - fr1 + 1;
if genOverlayImg
    %%
    hold on;
    ixSel = find((frmNo>=fr1) .* (frmNo<=fr2));
    ixSelTr = find((trInf2(:,1)>=fr1) .* ( (trInf2(:,1)+trInf2(:,2)) <=fr2));   % traces
    switch isTrace
        case 1 % spots
            xx = X(ixSel); yy = Y(ixSel);
            fr = frmNo(ixSel);
        case 2 % traces
            xx = TraceX(ixSel); yy = TraceY(ixSel);
            fr = frmNo(ixSel);
        case 3 % trace average
            xx = trInf2(ixSelTr,4); yy = trInf2(ixSelTr,5);
            fr = trInf2(ixSelTr,1);
        case 4 % trace recruit average (localized traces)
            ixSelTr2 = ixSelTr(trInf2(ixSelTr,6)>0.2);
            xx = trInf2(ixSelTr2,4); yy = trInf2(ixSelTr2,5);
            fr = trInf2(ixSelTr2,1);
    end

    hScat = scatter(xx,yy,10,CM(fr,:),'.');
    isPlay=1;
    if isPlay
        
        % figure
        figCtrl = figure;
        set(figCtrl,'Position',[10 100 450 200]);
        btx=0; bty=20;
        X= []; Y= X;
        hTextTime = uicontrol('style','text','BackgroundColor',[1 1 1],'String',sprintf('acqTime:%0.1fms, totTime:%.02fsec',acqTime*1000,acqTime*frames*binFrame),'Position',[btx+20 bty+150 250 15]);
        hTextSpotSel = uicontrol('style','text','BackgroundColor',[1 1 1],'String',sprintf('overlay data:%s',textDisp),'Position',[btx+20 bty+120 250 15]);
        hPopupData = uicontrol('Style', 'popup','BackgroundColor',[1 1 1],'Position', [btx+300 bty+150 150 25],...
            'String', {'spots','spots in traces','trace average','trace average (localized)'});    

        hFramesSlider1 = uicontrol('style','slider','units','pixel','position',[btx+280 bty+60 120 20]); 
        hFramesSlider2 = uicontrol('style','slider','units','pixel','position',[btx+280 bty+30 120 20]); 
        hNumFrames = uicontrol('style','slider','units','pixel','position',[btx+280 bty 120 20]); 
        hFramesSliderText1 = uicontrol('style','text','BackgroundColor',[1 1 1],'String',sprintf('fr1:%i',fr1),'units','pixel','position',[btx+200 bty+60 60 20]); 
        hFramesSliderText2 = uicontrol('style','text','BackgroundColor',[1 1 1],'String',sprintf('fr2:%i',fr2),'units','pixel','position',[btx+200 bty+30 60 20]); 
        hNumFramesSliderText = uicontrol('style','text','BackgroundColor',[1 1 1],'String',sprintf('#fr:%i',numFrames),'units','pixel','position',[btx+200 bty 60 20]); 

        hReset = uicontrol('style','pushbutton','BackgroundColor',[1 1 1],'String','reset','Position',[btx+20 bty+90 60 15],'Callback', 'updReset');
        hPrint = uicontrol('style','pushbutton','BackgroundColor',[1 1 1],'String','print','Position',[btx+100 bty+90 60 15],'Callback', 'callPrint');
        hPlay = uicontrol('style','pushbutton','BackgroundColor',[1 1 1],'String','play','Position',[btx+20 bty+60 160 15],'Callback', 'updPlay');
        hTextFrames = uicontrol('style','text','BackgroundColor',[1 1 1],'String','time:_(_-_)','Position',[btx+20 bty+30 160 15]);
        hPlaySpeedText2 = uicontrol('style','text','BackgroundColor',[1 1 1],'String','speed','units','pixel','position',[btx+20 bty 60 20]); 
        hplaySpeedSlider = uicontrol('style','slider','units','pixel','position',[btx+90 bty 90 20]); 

        %listen1 = addlistener(hFramesSlider1,'ActionEvent',@(hObject, event) updFrames(hObject, event,hScat,hPopupData,hTextFrames,hFramesSliderText1,hFramesSliderText2,hNumFramesSliderText,hFramesSlider1,hFramesSlider2,hNumFrames,hPlay,acqTime,X,Y,TraceX,TraceY,trInf2,frmNoSpot,frmNoTrace,CM,q,frames)); 
        listen1 = addlistener(hFramesSlider1,'ActionEvent',@(hObject, event) updFrames(hObject, event,hScat,hPopupData,hTextFrames,hFramesSliderText1,hFramesSliderText2,hNumFramesSliderText,hFramesSlider1,hFramesSlider2,hNumFrames,hPlay,acqTime,frames)); 
        %listen2 = addlistener(hFramesSlider2,'ActionEvent',@(hObject, event) updFrames(hObject, event,hScat,hPopupData,hTextFrames,hFramesSliderText1,hFramesSliderText2,hNumFramesSliderText,hFramesSlider1,hFramesSlider2,hNumFrames,hPlay,acqTime,X,Y,TraceX,TraceY,trInf2,frmNoSpot,frmNoTrace,CM,q,frames)); 
        listen2 = addlistener(hFramesSlider2,'ActionEvent',@(hObject, event) updFrames(hObject, event,hScat,hPopupData,hTextFrames,hFramesSliderText1,hFramesSliderText2,hNumFramesSliderText,hFramesSlider1,hFramesSlider2,hNumFrames,hPlay,acqTime,frames)); 
        %listen3 = addlistener(hNumFrames,'ActionEvent',@(hObject, event) updFrames(hObject, event,hScat,hPopupData,hTextFrames,hFramesSliderText1,hFramesSliderText2,hNumFramesSliderText,hFramesSlider1,hFramesSlider2,hNumFrames,hPlay,acqTime,X,Y,TraceX,TraceY,trInf2,frmNoSpot,frmNoTrace,CM,q,frames));
        listen3 = addlistener(hNumFrames,'ActionEvent',@(hObject, event) updFrames(hObject, event,hScat,hPopupData,hTextFrames,hFramesSliderText1,hFramesSliderText2,hNumFramesSliderText,hFramesSlider1,hFramesSlider2,hNumFrames,hPlay,acqTime,frames));
        listen4 = addlistener(hplaySpeedSlider,'ActionEvent',@(hObject, event) updPlaySpeed(hObject, event));
        listen5 = addlistener(hPopupData,'Value','PostSet',@(hObject, event) updFrames(hObject, event,hScat,hPopupData,hTextFrames,hFramesSliderText1,hFramesSliderText2,hNumFramesSliderText,hFramesSlider1,hFramesSlider2,hNumFrames,hPlay,acqTime,X,Y,TraceX,TraceY,trInf2,frmNoSpot,frmNoTrace,CM,q,frames)); 

        set(hFramesSlider1,'Min',1,'Max',frames);
        set(hFramesSlider2,'Min',1,'Max',frames);
        set(hNumFrames,'Min',1,'Max',frames);
        set(hplaySpeedSlider,'Min',0,'Max',frames);  % fps

        set(hFramesSlider1,'Value',fr1);
        set(hFramesSlider2,'Value',fr2);
        set(hNumFrames,'Value',fr1);

        % while loop
        q=0;
        framePlayTime=0;
        isPlaying = 0;
        numFrames = frames;
        %set(hNumFrames,'Value',fr1);
        while q == 0
            figure(figImg)
            btn = 0;      
            if ~isPlaying
                while btn == 0
                    btn = waitforbuttonpress;
                    k = get(figImg,'CurrentCharacter');
                end      
            end
            switch lower(k)
                case 'q'
                    q = 1;
            end

            fr1 = round(get(hFramesSlider1,'Value')); % current Frame
            fr2 = round(get(hFramesSlider1,'Value')); % last Frame
            %CM = genColorMap('jet',numFrames);
            if fr2 > frames, fr2 = frames; end

            % set data
            if isPlaying
                fr1 = df1+fr1;
                numFrames=dfnum+numFrames;
            end
            set(hFramesSlider1,'Value',fr1);
            set(hFramesSlider2,'Value',fr2);
            set(hNumFrames,'Value',fr1);

            t=round(numFrames*acqTime);
            t0=round(fr1*acqTime);
            t1=round(t0+t);
            set(hTextFrames,'String',sprintf('time:%i(%i-%i)  fr1=%i',t,t0,t1,fr1));

            % play speed
            playSpeed = get(hplaySpeedSlider,'Value');

            % frame process
            updScat;
            figure(figImg)
            colormap('gray')

            pause(framePlayTime);
        end

        %hTextImgMax = uicontrol('style','text','BackgroundColor',[1 1 1],'String','max','Position',[btx+20 bty+100 160 15]);
        %hImgMax = uicontrol('style','slider','units','pixel','position',[btx+20 bty+80 300 20]); 
        %hTextImgGamma = uicontrol('style','text','BackgroundColor',[1 1 1],'String','gamma','Position',[20 bty+140 160 15]);
        %hImgGamma = uicontrol('style','slider','units','pixel','position',[20 bty+120 300 20]); 

    end
    hold off;
elseif isDispRecruitment
    trSel = [];
    switch isTrace
        case 2 % traces
            xx = TraceX; yy = TraceY; int = TraceINT;
            fr = frmNoTrace;
            frInit = trInf(:,1); % starting frame
            if 0
                for i = 1:size(trInf,1)
                    trSel = [trSel trInf(i,3):trInf(i,3)+trInf(i,2)-1];
                end
                xx = xx(trSel);
                yy = yy(trSel);
                fr = fr(trSel);
                int = int(trSel);
            end
        case 3 % trace average (default)
            xx = trInf2(:,4); yy = trInf2(:,5);
            fr = trInf2(:,1);
        case 4 % trace recruit average (localized traces)
            ixSelTr2 = ixSelTr(trInf2(ixSelTr,6)>0.2);
            xx = trInf2(ixSelTr2,4); yy = trInf2(ixSelTr2,5);
            fr = trInf2(ixSelTr2,1);
    end    
    frmBin = floor(5/acqTime); % binning by 5 sec
%frmBin = 1;
    movieLen = frames*frmBinDetect*acqTime;
    frmLen = frmBin*frmBinDetect*acqTime; % [sec]
    
    % xyt convolution 
    if isTrace == 2
        tConvFrm = 30;
        xyConvPxHR = 1; % radius of the conv. ito highres. pix.
    elseif isTrace == 3 % trace average (default)
        tConvFrm = 6;
        xyConvPxHR = 1; % diameter of the conv. ito highres. pix.
    end
    tAvFrm = 6; % # of frames used in block average
    tAvSec = tAvFrm*frmLen; % [sec]
    tConvSec = tConvFrm*frmLen; % [sec]
    xyConvNm = xyConvPxHR*pxBin*pxSz;
    rateCoeff = 60/frmLen; % coeff forrecruitment per minute

    % round the spots on the edges
    imgSize = size(a);
    szY = imgSize(1);
    szX = imgSize(2);
    
    imgSizeBin = ceil(imgSize/pxBin);
    framesBin = ceil(frames/frmBin);
    binImg = zeros([imgSizeBin framesBin]); % binned image
    binImgRcrt = zeros([imgSizeBin framesBin]); % binned image for recruitment
    binImgRcrtIntW = zeros([imgSizeBin framesBin]); % binned image for recruitment
    binImgIntW  = zeros([imgSizeBin framesBin]); % intensity weighted
    binImgLong = binImg; % binned image keeping long traces
    binImgCodeByFrame = zeros(imgSizeBin); % binned image
    binImgRcrtTime = zeros(imgSizeBin); % binned image
    
    % color coding 
    isTimeCoding = 0;
    
    %% calculate intensity weighted trace coor.
    for i = 1:size(trInf,1)
        if (trInf(i,2) < minTrLen), continue; end;
        tracex = xx(trInf(i,3):trInf(i,3)+trInf(i,2)-1);
        tracey = yy(trInf(i,3):trInf(i,3)+trInf(i,2)-1);
        traceint = int(trInf(i,3):trInf(i,3)+trInf(i,2)-1);
        vs = ~isnan(tracex); % remove nans
        tracex = tracex(vs); tracey = tracey(vs); traceint = traceint(vs);
        trInf(i,9) = sum(tracex.*traceint)/sum(traceint); % weighted x position
        trInf(i,10) = sum(tracey.*traceint)/sum(traceint); % weighted y position
        trInf(i,11) = sum(traceint); % sum intensity
        xw = trInf(i,9); yw = trInf(i,9); 
        trInf(i,7) = (mean((tracex-xw).^2 + (tracey-yw).^2)).^0.5; % RMS displacement
        trInf(i,8) = findMaxDist(tracex,tracey); % MAX distance
    end
    
    %% delete traces (brightness and trace length and spread )
    mint = trInf(:,6); % mean intensity
    mint = sort(mint);
    nt = size(trInf,1); % number of traces
    intthresh = mint(round(nt*0.80)); % intensity threshold
    if ~intFilter, intthresh=0; end;
    nd = numel(xx); % data points
    ndelLowTrace = sum(trInf(:,6) < intthresh); % number of traces tb deleted
    ndelShortTrace = sum(trInf(:,2) < minTrLen); % number of traces tb deleted
    ndelSpreadTrace = sum(trInf(:,8) < cfg.trace.minXYspread); % number of traces tb deleted
    disp(sprintf('removing low intensity traces: %.02f%% of the traces are deleted',(ndelLowTrace)/size(trInf,1)*100 ));
    disp(sprintf('removing short traces: %.02f%% of the traces are deleted',(ndelShortTrace)/size(trInf,1)*100 ));
    disp(sprintf('removing spread traces: %.02f%% of the traces are deleted',(ndelSpreadTrace)/size(trInf,1)*100 ));
    i=1;
    hw = waitbar(0);
    while i <= size(trInf,1)
        if (trInf(i,2) < minTrLen) || (trInf(i,6) < intthresh) || trInf(i,8)<=cfg.trace.minXYspread 
            ixdel = [trInf(i,3)-1 trInf(i,3)+trInf(i,2)];
            trInf(i+1:end,3) = trInf(i+1:end,3) - trInf(i,2);
            if ixdel(1) == 0
                xx = xx(ixdel(2):end);
                yy = yy(ixdel(2):end);
                int = int(ixdel(2):end);
                fr = fr(ixdel(2):end);
                trInf = trInf(i+1:end,:);
            elseif ixdel(2) == numel(xx)
                xx = xx(1:ixdel(1));
                yy = yy(1:ixdel(1));
                int = int(1:ixdel(1));
                fr = fr(1:ixdel(1));
                trInf = trInf(1:i-1,:);
            else
                xx = [xx(1:ixdel(1)) xx(ixdel(2):end)];
                yy = [yy(1:ixdel(1)) yy(ixdel(2):end)];
                int = [int(1:ixdel(1)) int(ixdel(2):end)];
                fr = [fr(1:ixdel(1)) fr(ixdel(2):end)];
                trInf = [trInf(1:i-1,:); trInf(i+1:end,:)];
            end
        else
            i = i + 1;
        end
        waitbar(i/size(trInf,1),hw,'recruitmentTrack: filtering out traces...')
    end
    close(hw)
    disp(sprintf('removing low intensity traces: %.02f%% of the data points deleted',(nd-numel(xx))/nd*100 ));

    %trInf
        % 1: 1st frame
        % 2: number of frames
        % 3: position in the trace array
        % 4-6: mean x, y , int
        % 7 : RMS displacement
        % 8 : MAX distance
        % 9 : weighted x position
        % 10 : weighted y position
        % 11 : sum intensity
    
    % find(isnan(trInf(:,10))) % debug
%% BIN IMAGES    
    frTr = trInf(:,1); % frame of Traces
    offsetCM = ceil(framesBin/256); % offset for colormap
    hw = waitbar(0);
    for i = 1: framesBin % each binned frame
        fr1=(i-1)*frmBin+1; fr2=i*frmBin;
        ixSel = find((fr1<=fr) .* (fr<=fr2)); % selected coordnates
        trSel = find((fr1<=frTr) .* (frTr<=fr2)); % selected traces
        %ixSel2 = ; % selected traces
        %TraceXbin = xx(ixSel(xx(ixSel)>0));
        %TraceYbin = yy(ixSel(yy(ixSel)>0));
        TraceXbin = xx(ixSel);
        TraceYbin = yy(ixSel);
        TraceInt = int(ixSel);
        pxY = round(TraceYbin/pxBin);
        pxX = round(TraceXbin/pxBin);
        recXbin = trInf(trSel,9);
        recYbin = trInf(trSel,10);
        recInt = trInf(trSel,11);
        rX = round(recXbin/pxBin);
        rY = round(recYbin/pxBin);
        
        pxSel = (pxX>0) .* (pxX<=imgSizeBin(2)) .* (pxY>0) .* (pxY<=imgSizeBin(1));
        pxY = pxY(find(pxSel));
        pxX = pxX(find(pxSel));
        if isTimeCoding
            for j = 1:numel(pxX)
                binImgCodeByFrame(pxY(j),pxX(j))=i+offsetCM ...
                    +binImgCodeByFrame(pxY(j),pxX(j));
            end
        end
        if 0
            for j = 1:numel(pxX) % binning frames
                if 0 && trInf(j,2)>maxTrLen
                    binImgLong(pxY(j),pxX(j),i)=1 ...
                   +binImg(pxY(j),pxX(j),i);
                end
                binImg(pxY(j),pxX(j),i)=1 ...
               +binImg(pxY(j),pxX(j),i);
                binImgIntW(pxY(j),pxX(j),i)=TraceInt(j) ... % intensity weighted 
               +binImgIntW(pxY(j),pxX(j),i);
            end
        end
        for j = 1:numel(rY) % binning frames
            binImgRcrt(rY(j),rX(j),i)=1 ...
           +binImgRcrt(rY(j),rX(j),i);
            binImgRcrtIntW(rY(j),rX(j),i)=recInt(j) ...
           +binImgRcrtIntW(rY(j),rX(j),i);
            
        end       
        for j = 1:numel(rY) % binning frames
            binImgRcrtTime(rY(j),rX(j))=i+offsetCM ...
           +binImgRcrtTime(rY(j),rX(j));
        end       
        
        
        
        % intensity weighted trace coord.
        if 0 
            TraceXcoor = trInf(ixSel2,9);
            TraceYcoor = trInf(ixSel2,10);
            pxY = round(TraceYcoor/pxBin);
            pxX = round(TraceXcoor/pxBin);
            pxSel = (pxX>0) .* (pxX<=imgSizeBin(2)) .* (pxY>0) .* (pxY<=imgSizeBin(1));
            pxY = pxY(find(pxSel));
            pxX = pxX(find(pxSel));
        end
%         for j = 1:numel(pxX) % binning frames
%             binImgW(pxY(j),pxX(j),i)=1 ...
%            +binImgW(pxY(j),pxX(j),i);
%         end
        waitbar(i/framesBin,hw,'recruitmentTrack: binning images...')
    end
    close(hw)
    %binImg = binImgW; disp('using intensity weighted trace coord.');
    %binImgCodeByFrame = binImgCodeByFrame./sum(binImg,3);
    binImgRcrtTime = binImgRcrtTime./sum((binImgRcrt),3);
    
    %% xyt convolution
    funConv = fspecial('gaus',xyConvPxHR,xyConvPxHR/3);
    funConv(funConv<funConv(1,round(xyConvPxHR/2)))=0;
    funConv(funConv>0)=1;
    
    
    isCONVstack = 0;
    if isCONVstack 
        binImgConv = convn(binImg,repmat(funConv,[1 1 tConvFrm]),'full');
    end
    
    time = fix(clock);
    textTime = sprintf('%i/%i/%i, %ih%02im',time(1),time(2),time(3),time(4),time(5));
    
    textDataConv = 'recruitment per 5 sec';
    %% calc TIRF angle
    coeffTIRF_EPI_conversion = 0;
    TIRFstagePosEPI = cfg.acq.TIRFstagePosEPI;
    TIRFstagePosTIRF = cfg.acq.TIRFstagePosTIRF;
%    TIRFangle = calcTIRFangle(TIRFstagePosEPI,TIRFstagePosTIRF);
    TIRFangle = 0;
    %% calc bleaching power
    if 0
        bleachData = cfg.bleach.bleachData;
        nsPos = str2num(cell2mat(bleachData(:,1)));
        TIRFPos = str2num(cell2mat(bleachData(:,3)));
        for i = 1:size(bleachData,1)
            bleachTime_(i,:) = sscanf(cell2mat(bleachData(i,2)),'%im%is\n'); % minute and second
            bleachTime(i) = bleachTime_(i,1)*60+bleachTime_(i,2);
            bleachDataPrint{i,1} = nsPos(i);
            bleachDataPrint{i,2} = bleachTime(i);
            bleachDataPrint{i,3} = calcTIRFangle(TIRFstagePosEPI,TIRFPos(i));
        end
    end
    
    
    %% read cfg data
    if isfield(cfg.fit,'Coeff')
        Coeff               = cfg.fit.Coeff;
    else
        Coeff               = 0;
    end
    gausKernelSzPadded  = cfg.fit.detect.gausKernelSzPadded;
    gausKernelSg        = cfg.fit.detect.gausKernelSg;
    gausKernelSz        = cfg.fit.detect.gausKernelSz;
    gausHat1            = cfg.fit.detect.gausHat1;
    
    if isfield(cfg.fit,'gaus')
        tP2G_1  = cfg.fit.gaus.tP2G_1;
        tol     = cfg.fit.gaus.tol;
        sg      = cfg.fit.gaus.sg;
        sr      = cfg.fit.gaus.sr;
        intP    = cfg.fit.gaus.intP;
    end
    
    cellType    = cfg.cell.cellType;
    platingTime = cfg.cell.platingTime;
    surfTreat   = cfg.cell.surfTreat;
    chemTreat   = cfg.cell.chemTreat;
    cellOther   = cfg.cell.other;
    expDate     = cfg.exp.date;
    expTime     = cfg.exp.time;
    
    if isfield(cfg.trace,'traceJmpForCombination') && isTraceCombination
        traceJmpForCombination = cfg.trace.traceJmpForCombination;
    else traceJmpForCombination = 0; end;
    if isfield(cfg.trace,'sptJmpForTracing')
        sptJmpForTracing = cfg.trace.sptJmpForTracing;
    else sptJmpForTracing = cfg.trace.sptJmp; end;
    sptReAppearTime = cfg.trace.sptReAppearTime;
    maxTraceSpeed = cfg.trace.maxTraceSpeed;
    minTraceLength = cfg.trace.minTraceLength;   
    minXYspread = cfg.trace.minXYspread;

    tempBeforePositioning   = cfg.acq.tempBeforePositioning;
    tempChangeInTenMin      = cfg.acq.tempChangeInTenMin;
    power                   = cfg.acq.power;
    gain                    = cfg.acq.gain;
    NDfilt                  = cfg.acq.NDfilt; 
    
    if isTraceCombination
        isTraceCombinationYN = 'yes';
    else
        isTraceCombinationYN = 'no';
    end
    
    % coeffnorm
%     plot(CoeffNormAll); title('Coefficient change'); ylabel('Coefficient'); xlabel('time');
%     imgFig = getframe(gcf);
%     imgOut = imgFig.cdata;   
        
    
    %% generate TIFF info string
    tifDescriptBleach = 'NA';
    textProc1 = 'DETECTION';
    textProc2 = 'FIT-LOCALIZATION';
    textProc3 = 'TRACING';
    textProc4 = 'TRACE-SELECTION';
    textProc5 = 'PIXELATION';
    textProc5 = 'CONVOLUTION';
    textProc6 = 'BLOCK AVERAGE';
    tifDescriptProcess2 = 'Process2:NA';
    if 1% screen out
        tifDescriptImg = sprintf('\nIMG INFO: \nData\t\t\t: %s \nPixel Size\t\t: %.02fnm \nFrame Length\t: %.02fsec \nDuration\t\t: %.02fsec \nGeneration Time\t: %s \nfolder name\t\t: %s ',...
                                textDataConv,pxSz*pxBin,frmLen,movieLen,textTime,cd);
        tifDescriptExp = sprintf('\nEXP INFO: \nDate\t: %s \nTime\t: %s',expDate,expTime);
        tifDescriptCell = sprintf('\nCELL INFO: \nType\t\t\t\t: %s \nPlating Time\t\t: %s \nSurface Treatment\t: %s \nChemical Treatment\t: %s \nOther\t\t\t\t: %s', ...
                                cellType,platingTime,surfTreat,chemTreat,cellOther);
        tifDescriptAcq = sprintf('\nACQUISITION INFO: \nROI\t\t\t\t\t\t: %s \nPixel Size\t\t\t\t: %.02fnm \nacqTime\t\t\t\t\t: %.02fmsec \nTIRFangle\t\t\t\t: %.02f \ntempBeforePositioning\t: %s \ntempChangeIn10min\t\t: %s \nPower\t\t\t\t\t: %s \nGain\t\t\t\t\t: %s \nNDfilt\t\t\t\t\t: %s \nextraMagnification\t\t: %.01f',...
                                nmCrop,pxSz,acqTime*1000,TIRFangle,tempBeforePositioning,tempChangeInTenMin,power, gain, NDfilt, magEx);
        %tifDescriptBleach_ = sprintf('\nND filter position\t: %i \nbleach time\t\t\t: %isec \nTIRF angle\t\t\t: %.02f degrees',...
    %                            cell2mat(bleachDataPrint'));
        %tifDescriptBleach = sprintf('\nBLEACHING INFO:%s',  tifDescriptBleach_);
        tifDescriptProcess1 = sprintf('\nPROCESS 1\t\t: %s \nCoeff\t\t\t: %.02f \ngausKernelSg\t: %.02f pixels \nfrmBinDetect\t: %.02f', ...
                                                textProc1,Coeff,gausKernelSg,frmBinDetect);
        %tifDescriptProcess2 = sprintf('\nPROCESS 2\t\t\t: %s \nPosition Tolerance\t: %.02f pixels \nIntensity Tolerance\t: %.02f(multiplier) \nGaussian Sigma\t\t: %.02f \nGaussian Sigma Tolerance\t: %.02f (multiplier),\nPSF power(norm=1)\t: %.02f',...
        %                                        textProc2,tP2G_1,tol,sg,sr,intP);
        tifDescriptProcess3 = sprintf('\nPROCESS 3\t\t\t\t\t: %s \nMin. Trace Length\t\t\t: %i pixels \nSpot Re-appear Time\t\t\t: %i frames \nTrace Jump\t\t\t\t\t: %i pixels \nisTraceCombination\t\t\t: %s \nTrace Jump (Combination)\t: %i pixels \nTrace Selection (minTrLen)\t: %i frames (%.01fms) \nTrace Selection (xy spread)\t: %.02f px ', ...
                                                textProc3,minTraceLength,sptReAppearTime,sptJmpForTracing,isTraceCombinationYN,traceJmpForCombination,minTrLen,minTrLen*acqTime*1000,minXYspread);
        tifDescriptProcess4 = sprintf('\nPROCESS 4\t: %s \nfrmBin\t\t: %.02f \npxBin\t\t: %.02f',...
                                                textProc5,frmBin,pxBin);                                        
        tifDescriptProcess5 = sprintf('\nPROCESS 5\t: %s \nfrmBin\t\t: %.02f \npxBin\t\t: %.02f',...
                                                textProc5,frmBin,pxBin);
        tifDescriptProcess5 = sprintf('\nPROCESS 5\t: %s \nconvolution window (t)\t: %.02fsec (%i frames) \nconvolution window (xy)\t: %.02f nm (%i HR pixels) ',...
                                                textProc5,tConvSec,tConvFrm,xyConvNm,xyConvPxHR);
        tifDescriptProcess6 = sprintf('\nPROCESS 5\t: %s \nconvolution window (t)\t: %.02fsec (%i frames)',...
                                                textProc6,tAvSec,tAvFrm);                                        

        tifDescriptConv = sprintf('%s \n%s \n%s \n%s \n%s \n%s \n%s \n%s \n%s \n%s',tifDescriptImg,tifDescriptExp,tifDescriptCell,tifDescriptBleach,tifDescriptAcq,...
                    tifDescriptProcess1,tifDescriptProcess2,tifDescriptProcess3,tifDescriptProcess4,tifDescriptProcess5);
        tifDescriptBlockAv = sprintf('%s \n%s \n%s \n%s \n%s \n%s \n%s \n%s \n%s \n%s',tifDescriptImg,tifDescriptExp,tifDescriptCell,tifDescriptBleach,tifDescriptAcq,...
                    tifDescriptProcess1,tifDescriptProcess2,tifDescriptProcess3,tifDescriptProcess4,tifDescriptProcess6);            
        tifDescript = sprintf('%s \n%s \n%s \n%s \n%s \n%s \n%s \n%s \n%s \n%s',tifDescriptImg,tifDescriptExp,tifDescriptCell,tifDescriptBleach,tifDescriptAcq,...
                    tifDescriptProcess1,tifDescriptProcess2,tifDescriptProcess3,tifDescriptProcess4);
    end
    hfImgInf = fopen('imgInfo.txt','w');
    fprintf(hfImgInf,'%s \n%s \n%s \n%s \n%s \n%s \n%s \n%s \n%s \n%s',tifDescriptImg,tifDescriptExp,tifDescriptCell,tifDescriptBleach,tifDescriptAcq,...
                tifDescriptProcess1,tifDescriptProcess2,tifDescriptProcess3,tifDescriptProcess4,tifDescriptProcess5);
    fclose(hfImgInf);
    
    
    %% generate and write images
    
    % binImgColorSum
    %imwrite(uint16(binImgCodeByFrame),'binImgCodingTime.tif');
    %setTiffDescription('binImgCodingTime.tif',tifDescript);
    imwrite(uint16(binImgRcrtTime),'binImgRcrtTime.tif');
    setTiffDescription('binImgRcrtTime.tif',tifDescript);
    
    
    % Stacks and images to be generated 
    % 1 - binImg.tif : 5 second binning in time (sum)
    % 2 - binImg_blockSum30seconds.tif : 30 sec binning (sum)
    % 3 - binImg_blockSum30secondsColor.tif : 30 sec binning in color
    % 4 - binImgColorSum.tif : each recruitment coded by color for frame#
    % and summed
    % 5 - binImgConv.tif : 60 sec convolution of 5 second sums
    % 6 - binImgSum.tif : sum of all recruitment
    
    % binImg
    isBWstack = 0;
    scout = 0; % screen out
    if isBWstack
        binImg = uint16(binImg);
        fnBin = 'binImg.tif';
        stackWrite(binImg,fnBin);
        setTiffDescription(fnBin,tifDescript);
        if scout, tifDescript, end 
    end
    
    if isBWstack
        binImgIntW = uint16(binImgIntW/5);
        fnBin = 'binImgIntW.tif';
        stackWrite(binImgIntW,fnBin);
        setTiffDescription(fnBin,tifDescript);
        if scout, tifDescript, end 
    end
    
    
    if isBWstack
        binImgRcrt = uint16(binImgRcrt);
        fnBin = 'binImgRcrt.tif';
        stackWrite(binImgRcrt,fnBin);
        setTiffDescription(fnBin,tifDescript);
        if scout, tifDescript, end 
    end
    
    
    
    if isWrtLongTraces
        binImgLong = uint16(binImgLong);
        fnBinLong = 'binImgLong.tif';
        stackWrite(binImgLong,fnBinLong);
        setTiffDescription(fnBinLong,tifDescript);
    end
    
    
    
    % binImgConv
    if isCONVstack
        stackWrite(binImgConv,'binImgConv.tif');
        setTiffDescription('binImgConv.tif',tifDescriptConv)
        if scout, tifDescriptConv, end 
        
    end
    
    % binImgSum sum: BW image
    %binImgSum = sum(binImg,3);
    %imwrite(uint16(binImgSum),'binImgSum.tif');
    %setTiffDescription('binImgSum.tif',tifDescript)

    %binImgIntWSum = sum(double(binImgIntW)/1000,3);
    %imwrite(uint16(binImgIntWSum),'binImgIntWSum.tif');
    %setTiffDescription('binImgIntWSum.tif',tifDescript)

    binImgRcrtIntWSum = sum(double(binImgRcrtIntW)/1000,3);
    imwrite(uint16(binImgRcrtIntWSum),'binImgRcrtIntWSum.tif');
    setTiffDescription('binImgRcrtIntWSum.tif',tifDescript)

    binImgRcrtSum = sum(double(binImgRcrt),3);
    imwrite(uint16(binImgRcrtSum),'binImgRcrtSum.tif');
    setTiffDescription('binImgRcrtSum.tif',tifDescript)
    
    return;
    % binImg_blockSum30seconds block average
    nb = 5; % # frames in the block
    sy = size(binImg,1);
    sx = size(binImg,2);
    sz = size(binImg,3);
    % 1
    binImg2 = reshape(binImg,sy*sx,sz);
    binImg2 = reshape(binImg2(:,1:floor(sz/nb)*nb),sy*sx,nb,floor(sz/nb));
    binImg2 = sum(binImg2,2);
    binImgBlockAv = squeeze(binImg2);
    binImgBlockAv=reshape(binImgBlockAv,sy,sx,size(binImgBlockAv,2));
    % 2 : long
    binImg2 = reshape(binImgLong,sy*sx,sz);
    binImg2 = reshape(binImg2(:,1:floor(sz/nb)*nb),sy*sx,nb,floor(sz/nb));
    binImg2 = sum(binImg2,2);
    binImgLongBlockAv = squeeze(binImg2);
    binImgLongBlockAv = reshape(binImgLongBlockAv,sy,sx,size(binImgLongBlockAv,2));
    
    if size(binImgBlockAv,3)==0 
        display('not anough frames for binImg_blockSum30seconds');
        return;
    end
    stackWrite(binImgBlockAv,'binImg_blockSum30seconds.tif'); 
    if isWrtLongTraces && size(binImgLongBlockAv,3)>0, stackWrite(binImgLongBlockAv,'binImgLong_blockMax6frames.tif'); end;
    
    if 0 
        %% color the recruitment image
        fname_blockMax6frames = 'binImg_blockSum30secondsColor.tif';
        maxColor = 100;
        maxColor = max(binImgBlockAv(:));
        if max(binImgBlockAv(:)) > maxColor
            disp('ERROR: too much of recruitment')
        end

        binImgBlockAvColor = bw2color(binImgBlockAv,maxColor);

        for i = 1 : size(binImgBlockAv,3)
            if i == 1
                if exist(fname_blockMax6frames)
                    delete(fname_blockMax6frames);
                end
                imwrite(binImgBlockAvColor(:,:,:,i),fname_blockMax6frames,'Compression', 'none') 
            else
                imwrite(binImgBlockAvColor(:,:,:,i),fname_blockMax6frames,'WriteMode','append','Compression', 'none') 
            end
        end
        %max(max(binImgBlockAv,[],2),[],1)
        setTiffDescription('binImg_blockSum30seconds.tif',tifDescriptBlockAv,1);
        setTiffDescription('binImg_blockSum30secondsColor.tif',tifDescriptBlockAv,1);
    end
    return;
    %% debug color settings
    isTestColor = 0;
    if isTestColor
        N = 50;
        a = ones(N,N,3);
        a = a.*repmat(permute([1 0 1],[1 3 2]),N,N);
        a(1:20,1:20,:)=0;

        a4d = repmat(a,[1 1 1 5]);
        stackWriteColor(uint16(a4d.*2^16),'a4d.tif');
    end
else % 3D scatter
    XYT = [];
    for i = 1: frames
        ix=find(TraceX(:,i)>0);
        XYT = [XYT; [full(TraceX(ix,i)) full(TraceY(ix,i)) ix.*0+i]];
    end

end

if genImg 
    for i = 1:1 % frames     
        ix = find(xx(:,i)>0);
        if ~numel(ix), 
            %imwrite(zeros(size(IMG(:,:,1))),highResImg,'WriteMode','append');
            img = uint16(zeros(size(a)*resMag));
            if ~isOverSize, IMG(:,:,i) = img; end
            imwrite(img,highResImg,'WriteMode','append');
            continue; 
        end;

        frm = uint16(zeros(size(a)*resMag));
        pxX = round(xx(ix,i)*resMag); 
        pxY = round(yy(ix,i)*resMag);
        if isCountInc
            int_ = int(ix,i)>0;
        else    
            int_ = int(ix,i);
        end
        frm(sub2ind(size(frm),pxY,pxX)) = int_;
        img = reshape(frm,size(a)*resMag);
        img = logical(1-logical(img));
        if ~isOverSize, IMG(:,:,i) = img; end

        if i == 1
            delete(highResImg)
            imwrite(img,highResImg, 'Compression','packbits');
        else
            imwrite(img,highResImg,'WriteMode','append', 'Compression','packbits');
        end
        waitbar(i/frames)
    end
end

if genImg && ~isOverSize
    figure
    IMGproj = sum(IMG,3);
    imagesc(IMGproj);
end

%scatter(avX,avY,'.');

if genOverlayImg
    figure(figImg)
    labelTrLen=[];
    if isTrace>=2
        labelTrLen = sprintf('_minTrLen%.02fms',minTrLen*acqTime*1000);
    end
    labelTime = sprintf('%.02fms_%.02fsec_1stFrm%i',acqTime*1000,acqTime*numFrames,fr1);

    fOut = sprintf('recruitment%s%s_%s.tif',traceORspot{isTrace},labelTrLen,labelTime);
    imgOut = getframe(gcf);
    imwrite(imgOut.cdata,fOut,'Compression', 'none')
end

if exist('figCtrl'), close(figCtrl); end;
    save('cfg','cfg')
return
openfig('report.fig')
gfrm = getframe(gcf);
cfg.report.spots_imgSpy=gfrm.cdata;

%magImg = magImgPadded(R+1:end-R,R+1:end-R);
%imagesc(magImg);


%% debug


    
    TraceINT = TraceINT(1:6);
    TraceX = TraceX(1:6);
    TraceY = TraceY(1:6);
    frmNoTrace = frmNoTrace(1:6);
    trInf = trInf(1,:);
    save('traceData0-coeff-CMEanaly','TraceX','TraceY','TraceINT','frmNoTrace','trInf','cfg')
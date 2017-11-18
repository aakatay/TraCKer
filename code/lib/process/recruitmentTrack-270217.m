clear; close all;
F = findall(0,'type','figure','tag','TMWWaitbar'); delete(F)

acqTime = 0.03; %cfg.acq.acqTime;
% searches subfolders for trace data and display combined data
% run in cell folder with subfolders such as:
% tirf001_crop4\32X25Y44x40\_000-coeff762
% make sure there is only one data set per folder
fid = fopen('frames.txt','at');
fprintf(fid,'\nacqTime:%i',acqTime*1000);
fclose(fid);
 

isPrepOverlayData = 0; % debug: prepare using bin frame time=acq time

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
resMag = 1;
genImg = 0;
%if genImg, genOverlayImg = 0; end;
isCountInc = 1; % every recruitment counted as one
highResImg = 'highResImg.tif';
% pixelation
frmBinDetect = 1; % param used in tracker
pxBin = 0.25; % pixel size multiplier

fnXYZtraceData0 = rdir('**\_*\*traceData0-coeff*');
fnXYZ = rdir('**\_*\*xyzDataGaus-coeff*');
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

isComb = 1;
if isComb
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
else
    i=1;
    if isTrace>1, load(fnXYZtraceData0(i).name); end
    load(fnXYZ(i).name,'ixSptFrm');
end

% stats
if ~exist('.\stats'), mkdir('stats'), end;

pixelCAM                = 16; %str2num(cfg.acq.pixelCAM);
magObj                  = 100; %str2num(cfg.acq.magObj);
magEx                   = 1.5; %str2num(cfg.acq.magEx);

pxSz =  round(pixelCAM/magObj/magEx*10000)/10; %nm

if ~exist('cfg'), cfg = struct; end;
cfg.overlay = struct;
if ~isfield(cfg,'report'), cfg.report = struct; end;


% update minTrLen
minTrLen = 2;
maxTrLen = 10;

cfg.overlay.isCountInc = isCountInc;
%% data cropping 12X52Y47x47
isDataCropped=0; % = 1 if needs to be cropped
isNDcrop = 0; % ow imageJ crop
cfg.overlay.crop = struct;
cfg.overlay.crop.isNDcrop = isNDcrop;


fnm=rdir('**\*fname.mat');
load(fnm(1).name)


nCh = numel(fname);
last_= strfind(fname,'_');
if numel(last_) > 1
    nmCrop = fname(last_(end-2)+1:last_(end-1)-1);
    ixX = find(nmCrop=='X');
    ixY = find(nmCrop=='Y');
    ixx = find(nmCrop=='x');
    nmCrop2 = fname(last_(end-1)+1:last_(end)-1);
    ix_ = find(nmCrop2=='-');
end

if ~isempty(ixX) && ~isempty(ixY) && ~isempty(ixx)
    isDataCropped = 1;
    xCr = str2num(nmCrop(1:ixX-1));
    yCr = str2num(nmCrop(ixX+1:ixY-1));
    szXcr = str2num(nmCrop(ixY+1:ixx-1));
    szYcr = str2num(nmCrop(ixx+1:end));
    cfg.overlay.crop.xCr_X_yCr_Y_szXcr_x_szYcr = [xCr yCr szXcr szYcr]; 
    if ~isComb
        frame0 = str2num(nmCrop2(1:ix_-1));
        frames = str2num(nmCrop2(ix_+1:end))-frame0+1;
    end
end
cfg.overlay.crop.isDataCropped = isDataCropped;


maxNumFrames_time = 120; % 2min
if 120/acqTime > frames
    maxNumFrames_time = frames*acqTime;
end

if isempty(imgFile)
    a = zeros(szYcr,szXcr);
else
    a=imread(imgFile);
    if isDataCropped, a= a(yCr+1:yCr+szYcr,xCr+1:xCr+szXcr); end;
end

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

  
    cfg.report.disperseTraceElimRate = (size(trInf,1)-size(trInf2,1))/size(trInf,1);
end



switch isTrace
    case 1 % spots
        load(fnXYZ.name)
        textDisp = 'incidences';
        
        hist(INT,20);  title('intensity histogram'); pause;
        imgFig = getframe(gcf); imgOut = imgFig.cdata;
        imwrite(imgOut,'GausFit-intensity_histogram.tif');
        trInf = [];
        TraceINT = [];
    case 2 % traces
        textDisp = 'incidences in the traces';
    case 3 % trace average
        textDisp = 'trace average';
    case 4 % trace recruit average (localized traces)
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

trSel = [];
switch isTrace
    case 1
        xx = X; yy = Y; int = INT;
        fr = frmNoSpot;
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
binTime = 5; % binning by 5 sec
frmBin = floor(binTime/acqTime); 
frmBin2 = floor(binTime/10/acqTime); 
if isPrepOverlayData, frmBin = 1; end;
movieLen = frames*frmBinDetect*acqTime;
frmLen = frmBin*frmBinDetect*acqTime; % [sec]

% xyt convolution 
tConvFrm = 12;
xyConvPxHR = 1; % radius of the conv. ito highres. pix.

tAvFrm = 6; % # of frames used in block average
tAvSec = tAvFrm*frmLen; % [sec]
tConvSec = tConvFrm*frmLen; % [sec]
xyConvNm = xyConvPxHR*pxBin*pxSz;

% round the spots on the edges
imgSize = size(a);
szY = imgSize(1);
szX = imgSize(2);

imgSizeBin = ceil(imgSize/pxBin);
framesBin = ceil(frames/frmBin);
ThalfSec = 3*60; % [sec]
framesBin2 = ceil(ThalfSec/acqTime/frmBin2);
if isPrepOverlayData, framesBin = 1000; end
binImgRcrt = zeros([imgSizeBin framesBin]); % binned image for recruitment
binImgRcrtHalfSec = zeros([imgSizeBin framesBin]);
binImg = zeros([imgSizeBin framesBin]); % binned image
binImgRcrtIntW = zeros([imgSizeBin framesBin]); % binned image for recruitment
binImgRcrtTime = zeros(imgSizeBin); % binned image

sx = imgSizeBin(2);
sy = imgSizeBin(1);


%% calculate intensity weighted trace coor.
for i = 1:size(trInf,1)
    if (trInf(i,2) < minTrLen), continue; end;
    tracex = xx(trInf(i,3):trInf(i,3)+trInf(i,2)-1);
    tracey = yy(trInf(i,3):trInf(i,3)+trInf(i,2)-1);
    traceint = int(trInf(i,3):trInf(i,3)+trInf(i,2)-1);
    %vs = ~isnan(tracex); % remove nans
    vs = find(tracex>0);
    tracex = tracex(vs); tracey = tracey(vs); traceint = traceint(vs);
    mnx = mean(tracex); mny = mean(tracey); % mean position
    rx = tracex-mnx; ry = tracey-mny; % relative position
    trInf(i,9) = mnx+sum(rx.*traceint)/sum(traceint); % weighted x position
    trInf(i,10) = mny+sum(ry.*traceint)/sum(traceint); % weighted y position
    trInf(i,11) = sum(traceint); % sum intensity
    trInf(i,12) = max(traceint); % max intensity
    xw = trInf(i,9); yw = trInf(i,10); 
    trInf(i,7) = (mean((tracex-xw).^2 + (tracey-yw).^2)).^0.5; % RMS displacement
    if trInf(i,7)>5
        ccc=3; 
        
        %error('nananana')
        %break
    end
    trInf(i,8) = findMaxDist(tracex,tracey); % MAX displacement
end
%% display trace populations
%trInf
    % 1: 1st frame
    % 2: number of frames
    % 3: position in the trace array
    % 4-6: mean x, y , int
    % 7 : RMS displacement
    % 8 : MAX displacement
    % 9 : weighted x position
    % 10 : weighted y position
    % 11 : sum intensity
    % 12 : max intensity
%t11 = trInf(i,11);
if 1 
    figure(11)
    subplot(2,4,1);
    hist(trInf(:,12),100); title('max intensity')
    subplot(2,4,2);
    hist(trInf(:,6),100); title('mean intensity')
    subplot(2,4,3);
    hist(trInf(:,11),100); title('sum intensity')
    subplot(2,4,4);
    scatter(trInf(:,12),trInf(:,7),'.'); title('mean spread vs max intensity')
    subplot(2,4,5);
    scatter(trInf(:,12),trInf(:,11),'.'); title('sum intensity vs max intensity')
    subplot(2,4,6);
    scatter(trInf(:,12),trInf(:,6),'.'); title('mean intensity vs max intensity')
    subplot(2,4,7);
    scatter(trInf(:,8),trInf(:,7),'.'); title('mean spread vs max spread')
    subplot(2,4,8);
    scatter(trInf(:,8),trInf(:,2),'.'); title('trace length vs max spread')
    % write figure
    maximize
    imgFig = getframe(gcf);
    tpImg = imgFig.cdata; 
    tpFN = sprintf('trace populations.tif');
    imwrite(uint16(tpImg),tpFN)

    %% set2
    [N,X] = hist(trInf(:,1),100);

    figure(10)
    subplot(2,4,1);
    plot(X,smooth(N)); title('number of traces vs trace time')
    subplot(2,4,2);
    scatter(trInf(:,7),trInf(:,11),'.'); title('sum intensity vs mean spread')
    subplot(2,4,3);
    scatter(trInf(:,8),trInf(:,11),'.'); title('sum intensity vs max spread')
    subplot(2,4,4);
    scatter(trInf(:,1),trInf(:,12),'.'); title('max intensity vs trace time')
    subplot(2,4,5);
    scatter(trInf(:,1),trInf(:,6),'.'); title('mean intensity vs trace time')
    subplot(2,4,6);
    scatter(trInf(:,1),trInf(:,11),'.'); title('sum intensity vs trace time')
    subplot(2,4,7);
    scatter(trInf(:,8),trInf(:,7),'.'); title('mean spread vs max spread')
    subplot(2,4,8);
    scatter(trInf(:,8),trInf(:,2),'.'); title('trace length vs max spread')
    % write figure
    maximize
    imgFig = getframe(gcf);
    tpImg = imgFig.cdata; 
    tpFN = sprintf('trace populations2.tif');
    imwrite(uint16(tpImg),tpFN)
end
%return



%% delete traces (brightness and trace length and spread )

if isTrace > 1 % only if data is traces 
    %thrSpread = cfg.trace.minXYspread;
    thrSpread = 0.51; % spread threshold
    thrSpread = 10.51; % spread threshold
    %mxSpread = 0.8;
    selBrightestPercent = 100; % percentage of the selected traces acc. to brightness
    if selBrightestPercent==100, intFilter = 0;selBrightestPercent=90; else intFilter = 1;  end; % intensity filter
    delDim = (100-selBrightestPercent)/100;

    delShortTrace = find(trInf(:,2) < minTrLen); % number of traces tb deleted
    delLongTrace = find(trInf(:,2) > maxTrLen); % number of traces tb deleted
    ndelShortTrace = numel(delShortTrace);
    ndelLongTrace = numel(delLongTrace);
    mint = trInf(:,6); % mean intensity
    mint = sort(mint);
    nt = size(trInf,1); % number of traces
    intThresh = mint(round(nt*delDim)); % intensity threshold
    if ~intFilter, intThresh=0; end;
    nd = numel(xx); % data points
    delLowTrace = find(trInf(:,6) < intThresh); % number of traces tb deleted
    ndelLowTrace = numel(delLowTrace); % number of traces tb deleted
    delSpreadTrace = find(trInf(:,7) > thrSpread); % number of traces tb deleted
    ndelSpreadTrace = numel(delSpreadTrace);
    disp(sprintf('_________ long traces: %.02f%% of the traces are long and kept',(ndelLongTrace)/size(trInf,1)*100 ));
    disp(sprintf('removing short traces: %.02f%% of the traces are deleted',(ndelShortTrace)/size(trInf,1)*100 ));
    disp(sprintf('removing low intensity traces: %.02f%% of the traces are deleted',(ndelLowTrace)/size(trInf,1)*100 ));
    disp(sprintf('removing spread traces: %.02f%% of the traces are deleted',(ndelSpreadTrace)/size(trInf,1)*100 ));

    i=1;
    hw = waitbar(0);
    while i <= size(trInf,1)
        if ~isempty(find(i==delShortTrace)) || ~isempty(find(i==delLowTrace))|| ~isempty(find(i==delSpreadTrace))
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
            delShortTrace = delShortTrace-1;
            delLowTrace = delLowTrace-1;
            delSpreadTrace = delSpreadTrace-1;
        else
            i = i + 1;
        end
        waitbar(i/size(trInf,1),hw,'recruitmentTrack: filtering out traces...')
    end
    close(hw)
    disp(sprintf('removing traces: %.02f%% of the data points deleted',(nd-numel(xx))/nd*100 ));
else

    intThresh = 1700;
    ixNonFit = find(int<=intThresh);
    ixFit = find(int>intThresh);
    nonfitPercent = numel(ixNonFit)/numel(int)*100;
    disp(sprintf('GaussianFitLocalization: %.02f%% of the detections are removed',nonfitPercent))
    xx = xx(ixFit);
    yy = yy(ixFit);
    int = int(ixFit);
    fr = fr(ixFit);

end

%% display updated trace populations

if 0 
    close(11)

    figure(12)
    subplot(2,4,1);
    hist(trInf(:,12),100); title('max intensity')
    subplot(2,4,2);
    hist(trInf(:,6),100); title('mean intensity')
    subplot(2,4,3);
    hist(trInf(:,11),100); title('sum intensity')
    subplot(2,4,4);
    scatter(trInf(:,12),trInf(:,7),'.'); title('max intensity vs mean spread')
    subplot(2,4,5);
    scatter(trInf(:,12),trInf(:,11),'.'); title('max intensity vs sum intensity')
    subplot(2,4,6);
    scatter(trInf(:,12),trInf(:,6),'.'); title('max intensity vs mean intensity')
    subplot(2,4,7);
    scatter(trInf(:,8),trInf(:,7),'.'); title('max spread vs mean spread')
    subplot(2,4,8);
    scatter(trInf(:,8),trInf(:,2),'.'); title('max spread vs trace length')
    % write figure
    maximize
    imgFig = getframe(gcf);
    tpImg = imgFig.cdata; 
    tpFN = sprintf('trace populations-filtered.tif');
    imwrite(uint16(tpImg),tpFN)
end

% find(isnan(trInf(:,10))) % debug
%% BIN IMAGES    
if isTrace ~= 1, frTr = trInf(:,1); end; % frame of Traces
offsetCM = ceil(framesBin/256); % offset for colormap (1-frameShift)
hw = waitbar(0);
for i = 1: framesBin % each binned frame
    fr1=(i-1)*frmBin+1; fr2=i*frmBin; % 5 sec bins
    fr3=(i-1)*frmBin2+1; fr4=i*frmBin; % 0.5 sec bins
    ixSel = find((fr1<=fr) .* (fr<=fr2)); % selected coordnates
    TraceXbin = xx(ixSel);
    TraceYbin = yy(ixSel);
    TraceInt = int(ixSel);
    pxX = round(TraceXbin/pxBin);
    pxY = round(TraceYbin/pxBin);
    pxSel = (pxX>0) .* (pxX<=imgSizeBin(2)) .* (pxY>0) .* (pxY<=imgSizeBin(1));
    pxX = pxX(find(pxSel));    
    pxY = pxY(find(pxSel));
    TraceInt = TraceInt(find(pxSel));
    
    if isTrace > 1 % data is traces
        trSel = find((fr1<=frTr) .* (frTr<=fr2)); % selected traces (5sec bin)
        recXbin = trInf(trSel,9);
        recYbin = trInf(trSel,10);
        recInt = trInf(trSel,11);
        rX = round(recXbin/pxBin);
        rY = round(recYbin/pxBin);
        
        trSel = find((fr3<=frTr) .* (frTr<=fr4)); % selected traces (0.5sec bin)
        recXbin = trInf(trSel,9);
        recYbin = trInf(trSel,10);
        recInt = trInf(trSel,11);
        rX2 = round(recXbin/pxBin);
        rY2 = round(recYbin/pxBin);
    else
        rX = pxX;
        rY = pxY;
        recInt = TraceInt;
        pxX = [];
    end
        
    for j = 1:numel(rY) % binning frames
        if isnan(rX(j)) || rX(j)>sx || rY(j)>sy
            sprintf('exceeding localization = rX:%i, rY:%i',rX(j),rY(j));
            continue;
        end
        if i<=framesBin2
            binImgRcrtHalfSec(rY2(j),rX2(j),i)=1 ...
            +binImgRcrtHalfSec(rY2(j),rX2(j),i);
        end
   
        binImgRcrt(rY(j),rX(j),i)=1 ...
       +binImgRcrt(rY(j),rX(j),i);
        binImgRcrtIntW(rY(j),rX(j),i)=recInt(j) ...
       +binImgRcrtIntW(rY(j),rX(j),i);

        binImgRcrtTime(rY(j),rX(j))=i+offsetCM ...
       +binImgRcrtTime(rY(j),rX(j));
    end

      
    for j = 1:numel(pxX) % binning frames
        binImg(pxY(j),pxX(j),i)=1 ...
       +binImg(pxY(j),pxX(j),i);    
    end

    waitbar(i/framesBin,hw,'recruitmentTrack: binning images...')
end
close(hw)
binImgRcrtTime = binImgRcrtTime./sum((binImgRcrt),3);




time = fix(clock);
textTime = sprintf('%i/%i/%i, %ih%02im',time(1),time(2),time(3),time(4),time(5));

textDataConv = 'recruitment per 5 sec';


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
textProc5 = 'CONVOLUTION';
textProc6 = 'BLOCK AVERAGE';
tifDescriptProcess2 = 'Process2:NA';
TIRFangle = 0;
if 1% screen out
    tifDescriptImg = sprintf('\nIMG INFO: \nData\t\t\t: %s \nPixel Size\t\t: %.02fnm \nFrame Length\t: %.02fsec \nDuration\t\t: %.02fsec \nGeneration Time\t: %s \nfolder name\t\t: %s ',...
                            textDataConv,pxSz*pxBin,frmLen,movieLen,textTime,cd);
    tifDescriptExp = []; %sprintf('\nEXP INFO: \nDate\t: %s \nTime\t: %s',expDate,expTime);
    tifDescriptCell = []; %sprintf('\nCELL INFO: \nType\t\t\t\t: %s \nPlating Time\t\t: %s \nSurface Treatment\t: %s \nChemical Treatment\t: %s \nOther\t\t\t\t: %s', ...
                           %cellType,platingTime,surfTreat,chemTreat,cellOther);
    tifDescriptAcq = []; %sprintf('\nACQUISITION INFO: \nROI\t\t\t\t\t\t: %s \nPixel Size\t\t\t\t: %.02fnm \nacqTime\t\t\t\t\t: %.02fmsec \nTIRFangle\t\t\t\t: %.02f \ntempBeforePositioning\t: %s \ntempChangeIn10min\t\t: %s \nPower\t\t\t\t\t: %s \nGain\t\t\t\t\t: %s \nNDfilt\t\t\t\t\t: %s \nextraMagnification\t\t: %.01f',...
                            %nmCrop,pxSz,acqTime*1000,TIRFangle,tempBeforePositioning,tempChangeInTenMin,power, gain, NDfilt, magEx);
    tifDescriptProcess1 = []; %sprintf('\nPROCESS 1\t\t: %s \nCoeff\t\t\t: %.02f \ngausKernelSg\t: %.02f pixels \nfrmBinDetect\t: %.02f', ...
                             %               textProc1,Coeff,gausKernelSg,frmBinDetect);
    %tifDescriptProcess2 = sprintf('\nPROCESS 2\t\t\t: %s \nPosition Tolerance\t: %.02f pixels \nIntensity Tolerance\t: %.02f(multiplier) \nGaussian Sigma\t\t: %.02f \nGaussian Sigma Tolerance\t: %.02f (multiplier),\nPSF power(norm=1)\t: %.02f',...
    %                                        textProc2,tP2G_1,tol,sg,sr,intP);
    tifDescriptProcess3 = []; %sprintf('\nPROCESS 3\t\t\t\t\t: %s \nMin. Trace Length\t\t\t: %i pixels \nSpot Re-appear Time\t\t\t: %i frames \nTrace Jump\t\t\t\t\t: %i pixels \nisTraceCombination\t\t\t: %s \nTrace Jump (Combination)\t: %i pixels \nTrace Selection (minTrLen)\t: %i frames (%.01fms) \nTrace Selection (xy spread)\t: %.02f px ', ...
                               %             textProc3,minTraceLength,sptReAppearTime,sptJmpForTracing,isTraceCombinationYN,traceJmpForCombination,minTrLen,minTrLen*acqTime*1000,minXYspread);
    tifDescriptProcess4 = sprintf('\nPROCESS 4\t: %s \nfrmBin\t\t: %.02f \npxBin\t\t: %.02f',...
                                            textProc5,frmBin,pxBin);                                        
    tifDescriptProcess5 = sprintf('\nPROCESS 5\t: %s \nconvolution window (t)\t: %.02fsec (%i frames) \nconvolution window (xy)\t: %.02f nm (%i HR pixels) ',...
                                            textProc5,tConvSec,tConvFrm,xyConvNm,xyConvPxHR);
    tifDescriptProcess6 = sprintf('\nPROCESS 6\t: %s \nconvolution window (t)\t: %.02fsec (%i frames)',...
                                            textProc6,tAvSec,tAvFrm);                                        
                                        
    tifDescriptConv = sprintf('%s \n%s \n%s \n%s \n%s \n%s \n%s \n%s \n%s \n%s',tifDescriptImg,tifDescriptExp,tifDescriptCell,tifDescriptBleach,tifDescriptAcq,...
                tifDescriptProcess1,tifDescriptProcess2,tifDescriptProcess3,tifDescriptProcess4,tifDescriptProcess5);
    tifDescriptBlockAv = sprintf('%s \n%s \n%s \n%s \n%s \n%s \n%s \n%s \n%s \n%s',tifDescriptImg,tifDescriptExp,tifDescriptCell,tifDescriptBleach,tifDescriptAcq,...
                tifDescriptProcess1,tifDescriptProcess2,tifDescriptProcess3,tifDescriptProcess4,tifDescriptProcess6);            

    tifDescript = sprintf('%s \n%s \n%s \n%s \n%s \n%s \n%s \n%s \n%s ',tifDescriptImg,tifDescriptExp,tifDescriptCell,tifDescriptBleach,tifDescriptAcq,...
                tifDescriptProcess1,tifDescriptProcess2,tifDescriptProcess3);
end
hfImgInf = fopen('imgInfo.txt','w');
fprintf(hfImgInf,'%s \n%s \n%s \n%s \n%s \n%s \n%s \n%s \n%s \n%s',tifDescriptImg,tifDescriptExp,tifDescriptCell,tifDescriptBleach,tifDescriptAcq,...
            tifDescriptProcess1,tifDescriptProcess2,tifDescriptProcess3);
fclose(hfImgInf);


%% generate and write images

imwrite(uint16(binImgRcrtTime==9),'binImgRcrtTime.tif');
imwrite(uint16(binImgRcrtTime),'binImgRcrtTime.tif');
setTiffDescription('binImgRcrtTime.tif',tifDescript);


%imwrite(uint16(binImgRcrtTime(4*yCr+1:4*yCr+4*szY,4*xCr+1:4*xCr+4*szX)==9),'binImgRcrtTimeBUG.tif');
% Stacks and images to be generated 

% 
scout = 0; % screen out
sz = size(binImg,3);
%if ~isPrepOverlayData
    binImgRcrt = uint16(binImgRcrt);
    fnBin = 'binImgRcrt.tif';
    stackWrite(binImgRcrt,fnBin);
    setTiffDescription(fnBin,tifDescript);
    if scout, tifDescript, end 
    
    
    binImgRcrtHalfSec = uint16(binImgRcrtHalfSec);
    fnBin = 'binImgRcrtHalfSec.tif';
    stackWrite(binImgRcrtHalfSec,fnBin);
    setTiffDescription(fnBin,tifDescript);
    if scout, tifDescript, end 

    %% images of recruitment

    binImgRcrtIntWSum = sum(double(binImgRcrtIntW)/1000,3);
    imwrite(uint16(binImgRcrtIntWSum),'binImgRcrtIntWSum.tif');
    setTiffDescription('binImgRcrtIntWSum.tif',tifDescript);

    binImgRcrtSum = sum(double(binImgRcrt),3);
    imwrite(uint16(binImgRcrtSum),'binImgRcrtSum.tif');
    setTiffDescription('binImgRcrtSum.tif',tifDescript);
    recruitmentRate = sum(reshape(double(binImgRcrt),sx*sy,sz),1)/5; % per sec
    %clear binImgRcrt;
%end
%return;
%% convolution/blockSum stacks of detections
    if isPrepOverlayData, return; end
    display('convolution/blockSum stacks of detections')
    binImg = uint16(binImg);
    fnBin = 'binImg.tif';
    stackWrite(binImg,fnBin);
    % binImgConv
    %% xyt convolution
    funConv = fspecial('gaus',xyConvPxHR,xyConvPxHR/3);
    funConv(funConv<funConv(1,round(xyConvPxHR/2)))=0;
    funConv(funConv>0)=1;
    
    isCONVstack = 1;
    if isCONVstack 
        binImgConv1min = convn(binImg,repmat(funConv,[1 1 12]),'full');
        stackWrite(binImgConv1min,'binImgConv1min.tif');
        setTiffDescription('binImgConv1min.tif',tifDescriptConv)
        if scout, tifDescriptConv, end 
        clear binImgConv1min
        
%         binImgRcrtConv1min = convn(binImgRcrt,repmat(funConv,[1 1 12]),'full');
%         stackWrite(binImgRcrtConv1min,'binImgRcrtConv1min.tif');
%         setTiffDescription('binImgRcrtConv1min.tif',tifDescriptConv)
%         if scout, tifDescriptConv, end 
%         clear binImgRcrtConv1min

        binImgRcrtConv1min = convn(binImgRcrt,repmat(funConv,[1 1 12]),'full');
        stackWrite(binImgRcrtConv1min,'binImgRcrtConv1min.tif');
        setTiffDescription('binImgRcrtConv1min.tif',tifDescriptConv)
        if scout, tifDescriptConv, end 
        clear binImgRcrtConv1min
        
        binImgRcrtConv3sec = convn(binImgRcrtHalfSec,repmat(funConv,[1 1 6]),'full');
        stackWrite(binImgRcrtConv3sec,'binImgRcrtConv3sec.tif');
        setTiffDescription('binImgRcrtConv3sec.tif',tifDescriptConv)
        if scout, tifDescriptConv, end 
        clear binImgRcrtConv3sec
    end
    
    % binImg_blockSum30seconds block average
    nb = tAvFrm; % # frames in the block
    % 1
    binImg2 = reshape(binImg,sy*sx,sz);
    binImg2 = reshape(binImg2(:,1:floor(sz/nb)*nb),sy*sx,nb,floor(sz/nb));
    binImg2 = sum(binImg2,2);
    binImgBlockAv = squeeze(binImg2);
    binImgBlockAv=reshape(binImgBlockAv,sy,sx,size(binImgBlockAv,2));
    stackWrite(binImgBlockAv,'binImg_blockSum30seconds.tif'); 
    setTiffDescription('binImg_blockSum30seconds.tif',tifDescriptBlockAv,1);
    clear binImgBlockAv
%% photobleaching
close all;
    % calculate half bleach time
    mxrecRate=max(recruitmentRate);
    recruitmentRateSmth = smooth(recruitmentRate);
    hlfRecRate=mxrecRate/2;
    ix1 = find(recruitmentRateSmth< hlfRecRate,1);
    v1=recruitmentRateSmth(ix1-1);
    v2=recruitmentRateSmth(ix1);
    hlfBlchTime = (ix1+(v2-hlfRecRate)/(v1-v2))*binTime;
%plot
figure(13)
xTime = (1:sz)*binTime; % [sec]
plot(xTime,recruitmentRate); 
ax = gca;
ax.XTick = [(0:3:sz)*binTime];
grid minor;
maximize
ylabel('# of recruitment incidences/sec')
xlabel('time [sec]');
title(sprintf('recruitmentRate[incidence/sec]. half bleach time: %.02fsec',hlfBlchTime))
imgFig = getframe(gcf);
recRateImg = imgFig.cdata; 
rrFN = sprintf('recruitmentRate-hlfTime%.02fsec.tif',hlfBlchTime);
imwrite(uint16(recRateImg),rrFN)
pause
close(13)

return;


%% debug


    
    TraceINT = TraceINT(1:6);
    TraceX = TraceX(1:6);
    TraceY = TraceY(1:6);
    frmNoTrace = frmNoTrace(1:6);
    trInf = trInf(1,:);
    save('traceData0-coeff-CMEanaly','TraceX','TraceY','TraceINT','frmNoTrace','trInf','cfg')
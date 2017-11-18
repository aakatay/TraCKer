clear; close all;
PWD = pwd;
ix=strfind(PWD,'\');
PWD = PWD(ix(end)+1:end);
if strcmp(PWD(1),'_'), cd('..'); end;
F = findall(0,'type','figure','tag','TMWWaitbar'); delete(F)

acqTime = 0.03; %cfg.acq.acqTime;
% searches subfolders for trace data and display combined data
% run in cell folder with subfolders such as:
% tirf001_crop4\32X25Y44x40\_000-coeff762
% make sure there is only one data set per folder
fid = fopen('frames.txt','at');
fprintf(fid,'\nacqTime:%i',acqTime*1000);
fclose(fid);

isTrace=2; % ============================ select overlay data =====================================
%1:spots 2:traces 3:trace average 4:trace recruit average (localized traces) 

%% load files
% parameters
highResImg = 'highResImg.tif';
% pixelation
frmBinDetect = 1; % param used in tracker
pxBin = 0.25; % pixel size multiplier

fnXYZtraceData0 = rdir('**\_*\*traceData0-coeff*');
fnXYZ = rdir('**\_*\*xyzDataGaus-coeff*');

TraceX_ = [];
TraceY_ = [];
TraceINT_ = [];
frmNoTrace_ = [];
trInf_ = [];
frames = 0;
tracePos = 0;

load(fnXYZtraceData0(1).name);
load(fnXYZ(1).name,'ixSptFrm');

% stats
if ~exist('.\stats'), mkdir('stats'), end;

% update minTrLen
minTrLen = 2;
maxTrLen = 10000;

%% data cropping 12X52Y47x47
isDataCropped=0; % = 1 if needs to be cropped
isNDcrop = 0; % ow imageJ crop

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
    frame0 = str2num(nmCrop2(1:ix_-1));
    frames = str2num(nmCrop2(ix_+1:end))-frame0+1;
end
a = zeros(szYcr,szXcr);
%% trace - spot loop
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

% correct TraceInt values
TraceInt = TraceINT; % ---
TraceInt(TraceInt<0)=0; TraceInt(isnan(TraceInt))=0;


fr1 = 1;
fr2 = frames;
numFrames = fr2 - fr1 + 1;

trSel = [];
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

binTime = 5; % binning by 5 sec
frmBin = (binTime/acqTime);  

% xyt convolution 
xyConvPxHR = 1; % radius of the conv. ito highres. pix.


% round the spots on the edges
imgSize = size(a);
szY = imgSize(1);
szX = imgSize(2);

imgSizeBin = ceil(imgSize/pxBin);
framesBin = ceil(frames/frmBin);
binImgRcrt = zeros([imgSizeBin framesBin]); % binned image for recruitment
binImgRcrtHalfSec = zeros([imgSizeBin framesBin]);
binImg = zeros([imgSizeBin framesBin]); % binned image
binImgRcrtIntW = zeros([imgSizeBin framesBin]); % binned image for recruitment
binImgRcrtTime = zeros(imgSizeBin); % binned image

sx = imgSizeBin(2);
sy = imgSizeBin(1);


%% calculate intensity weighted trace coor.
for i = 1:size(trInf,1)
    if (trInf(i,2) < minTrLen), continue; end
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
    clear xw yw;
    trInf(i,8) = findMaxDist(tracex,tracey); % MAX displacement
end

FN = 'traceData_recTrack.mat';
save(FN,'acqTime','TraceX', 'TraceY', 'TraceINT', 'frmNoTrace', 'trInf','ixSptFrm','cfg');

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
recruitmentTrackDisp1; % display trace stats

%thrSpread = 1; % spread threshold [px];
%recruitmentTrack_FilterTraces; 

%% BIN IMAGES    
if isTrace ~= 1, frTr = trInf(:,1); end; % frame of Traces
offsetCM = ceil(framesBin/256); % offset for colormap (1-frameShift)
hw = waitbar(0);
for i = 1: framesBin % each binned frame
    %fr1=(i-1)*frmBin+1; fr2=i*frmBin; % 5 sec bins
    fr1=floor((i-1)*frmBin)+1; fr2=floor(i*frmBin); % 5 sec bins 10/19/2017	
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
    
    trSel = find((fr1<=frTr) .* (frTr<=fr2)); % selected traces (5sec bin)
    recXbin = trInf(trSel,9);
    recYbin = trInf(trSel,10);
    recInt = trInf(trSel,11);
    rX = round(recXbin/pxBin);
    rY = round(recYbin/pxBin);
        
    for j = 1:numel(rY) % binning frames
        if isnan(rX(j)) || rX(j)>sx || rY(j)>sy
            sprintf('exceeding localization = rX:%i, rY:%i',rX(j),rY(j));
            continue;
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

%% generate and write images

imwrite(uint16(binImgRcrtTime==9),'binImgRcrtTime.tif');
imwrite(uint16(binImgRcrtTime),'binImgRcrtTime.tif');

% 
sz = size(binImg,3);
binImgRcrt = uint16(binImgRcrt);
fnBin = 'binImgRcrt.tif';
stackWrite(binImgRcrt,fnBin);



%% images of recruitment
binImgRcrtIntWSum = sum(double(binImgRcrtIntW)/1000,3);
imwrite(uint16(binImgRcrtIntWSum),'binImgRcrtIntWSum.tif');

binImgRcrtSum = sum(double(binImgRcrt),3);
imwrite(uint16(binImgRcrtSum),'binImgRcrtSum.tif');
recruitmentRate = sum(reshape(double(binImgRcrt),sx*sy,sz),1)/5; % per sec

%% time convolution
binImgRcrtConv1min = convn(binImgRcrt,ones(1,1,12),'full');
stackWrite(binImgRcrtConv1min,'binImgRcrtConv1min.tif');

%% convolution/blockSum stacks of detections
display('convolution/blockSum stacks of detections')
binImg = uint16(binImg);
fnBin = 'binImg.tif';
stackWrite(binImg,fnBin);

    
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
close(13)

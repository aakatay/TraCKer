% displays the localization and traces in an overlay movie
% generates the movie by concatenates crops in a grid 
% uses CMEpitFinder info for crops

% E:\MATLAB\SUM-AP2\170605-SUM-AP2\cell14
clear all;
close all;
F = findall(0,'type','figure','tag','TMWWaitbar'); delete(F)
acqTime = 10e-3; %acquisition time of acq movie
lastFrm0 = 10;
%lastFrm0 = 0;
% 
windowSz = 20; % crop window size for display
windowSz0 = 6; % smaller for only pit data
szpad2 = 2; % pad to separate crops
sz = szpad2*2+windowSz; % window size with pads

%% files
acqImgFN='acqPitCrops.tif';
load pitCoors; % XY LTpits

if 0
    movefile('pitCoors.mat','pitCoors0.mat');
    pitSel = [1 2 3 9 18 21];
    XY = XY(pitSel,:);
    save('pitCoors','XY');
end
traceMAT = rdir('**\_*\traceData0*.mat');
traceMAT = traceMAT(1).name;
load(traceMAT);

% output file
pitRecruitmentProfileFN = 'pitRecruitmentProfile2.tif'; % rec time profiles

%% pit info
np = size(XY,1); % number of pits
XY = round(XY/4);
x = XY(:,1);
y = XY(:,2);
w = windowSz/2;
w0 = windowSz0/2;

rX = trInf(:,4); % recruitment data
rY = trInf(:,5);

%% loop to update to crop data
strInt = [];
trInf2 = [];
fr = zeros(1e5,1);
for i = 1:np % each pit
    Lt = LtPits(:,:,i);
    ixs = find(Lt(sub2ind(size(Lt),round(rY*4),round(rX*4))));
    
    fr = trInf(ixs,1);
    FR(fr,i) = 1; % recruitment times (binary)
    %Int(i) = numel(ixs0); % intensity of the pit
    trInf2_ = trInf(ixs,:);
    trInf2_(:,9) = - x(i)+w+1; % x-offset
    trInf2_(:,10) = - y(i)+w+1; % y-offset
    trInf2_(:,4) = trInf2_(:,4) - x(i)+w+1 + (i-1)*sz;
    trInf2_(:,5) = trInf2_(:,5) - y(i)+w+1;
    trInf2_(:,1) = trInf2_(:,1) - fr(1)+1; % frst frame
    
    frstFrame(i) = fr(1);
    trInf2 = [trInf2; trInf2_];
    %strInt = sprintf('%s,''No%i,int:%i''',strInt,i,Int(i));
end

%imagesc(conv2(FR',ones(1,20)))
if lastFrm0 == 0
    lastFrm = find(sum(FR,2)>0,1,'last');
else
    lastFrm = lastFrm0;
end

FR = FR(1:lastFrm,:);
FRshft = repmat([1:np]-1,lastFrm,1);
FR = FR + FRshft;

% save 
trInf = trInf2;
save('traceDataCrops','trInf')

%% edit the acquisition movie
save('pitCropsInfo','frstFrame','lastFrm','windowSz','szpad2'); % read by combineROI
combineROI; % generating acqPitCrops.tif

%% 
iminf = imfinfo('acqPitCrops.tif');
nf =numel(iminf); % number of frames
dispTracesCore(nf,1);

%% plot trace time profiles
if 0 
    figure(1000)
    maximize;
    t = (1:lastFrm)*acqTime;
    plot(t,FR,'Linewidth',2)
    strInt(1) = [];
    legendCall = sprintf('legend(%s)',strInt);
    eval(legendCall)

    title({'recruitment number vs time';'legend names are total number of recruitments for each pit'})
    xlabel('time (sec)')
    ylabel('# recruitments')


    imgFig = getframe(gcf); 
    figCap_recProf = imgFig.cdata;
    imwrite(figCap_recProf,pitRecruitmentProfileFN,'Compression', 'none') 
end
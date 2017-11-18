% displays the localization and traces in an overlay movie
% generates the movie by concatenates crops in a grid 
% uses CMEpitFinder info for crops
% generates 'acqPitCrops.tif' & 'TraceImage-acqPitCrops.tif'

% E:\MATLAB\SUM-AP2\170605-SUM-AP2\cell14
clear all;
close all;
F = findall(0,'type','figure','tag','TMWWaitbar'); delete(F)
mxMovieLen = 1000;
mxMovieLen = 0;
%lastFrm0 = 0;
% 
windowSz = 20; % crop window size for display
windowSz0 = 6; % smaller for only pit data
szpad2 = 2; % pad to separate crops
sz = szpad2*2+windowSz; % window size with pads

%% files
load pitCoors; % XY LTpits

if 0
    movefile('pitCoors.mat','pitCoors0.mat');
    pitSel = [1 2 3 9 18 21];
    XY = XY(pitSel,:);
    save('pitCoors','XY');
end
% output file
pitRecruitmentProfileFN = 'pitRecruitmentProfile2.tif'; % rec time profiles

%% data
% TraceX TraceY TraceX2 TraceY2 
fn2 = rdir('**\traceData0-coeff*.mat');
fn3 = rdir('**\traceJmplessData-coeff*.mat');
load(fn2(1).name);
load(fn3(1).name);

% trInf
load('traceData_recTrack.mat');

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
for i = 1:np % each pit
    Lt = LtPits(:,:,i);
    LtWide = conv2(Lt,ones(3),'same');
    ixs = find(Lt(sub2ind(size(Lt),round(rY*4),round(rX*4))));
    ixs = find(LtWide(sub2ind(size(Lt),round(rY*4),round(rX*4))));
    
    fr = trInf(ixs,1);
    FR(fr,i) = 1; % recruitment times (binary)
    %Int(i) = numel(ixs0); % intensity of the pit
    trInf2_ = trInf(ixs,:);
    trInf2_(:,13) = - x(i)+w+1 + (i-1)*sz + szpad2; % x-offset
    trInf2_(:,14) = - y(i)+w+1 + szpad2; % y-offset
    trInf2_(:,9) = trInf2_(:,9) + trInf2_(:,13);
    trInf2_(:,10) = trInf2_(:,10) + trInf2_(:,14);
    trInf2_(:,1) = trInf2_(:,1) - fr(1)+1; % frst frame
    
    frstFrame(i) = fr(1);
    trInf2 = [trInf2; trInf2_];
    %strInt = sprintf('%s,''No%i,int:%i''',strInt,i,Int(i));
    lastFrm(i) = find(FR(:,i)>0,1,'last');
end

%imagesc(conv2(FR',ones(1,20)))
if mxMovieLen == 0
    mxMovieLen = max(lastFrm-frstFrame);
end
% save 
trInf = trInf2;
save('traceDataCrops','trInf')

%% edit the acquisition movie
save('pitCropsInfo','frstFrame','lastFrm','mxMovieLen','windowSz','szpad2'); % read by combineROI
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
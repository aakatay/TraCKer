% reads and overlays two channels
% run in lower folder
%ch1: 3x3 center empty
%ch2: 3x3 only center
clear

load fname;
fname='acqDiff_0X0Y512x512_1-255_006.tif';
%fname='acqDiff_0X0Y512x512_1-255_006_235X46Y111x96.tif';

%fname='acqDiff_0X0Y512x512_1-354_003.tif';
fnameOrig=fname;
fnameOrig(4:7)=[]; % remove 'Diff'

shft = 0;
shft1 = 2;
shft2 = 3;

fnameCLCorig = sprintf('clc\\%s',fnameOrig);
fnameCLC = sprintf('clc\\%s',fname);
fnameDYN = sprintf('dyn\\%s',fname);
overlayFN = 'ClcDyn3x3Overlay.tif';

infFN=imfinfo(fnameDYN);
Frames = numel(infFN);
infFN=imfinfo(fnameCLC);
Frames = numel(infFN);

A_ = imread(fnameCLC,1); % clathrin channel
mag=3;
sz=size(A_);
szx=sz(2)*mag;
szy=sz(1)*mag;
%Frames=10;

if exist(overlayFN), delete(overlayFN); end;
hw =  waitbar(0);
for i = 1:Frames
    A_ = imread(fnameCLC,i); % clathrin channel
    A_ = repelem(A_,mag,mag);
    D_ = imread(fnameDYN,i+shft1); % dynamin channel
    D_ = repelem(D_,mag,mag);
    Aorig = imread(fnameCLCorig,i+shft2);
    Aorig = repelem(Aorig,mag,mag);
    
    % only in the center
    D = D_;
    D2_ = zeros(size(D_));
    D(:,1:3:szx)=0;
    D(:,3:3:szx)=0;
    D(1:3:szy,:)=0;
    D(3:3:szy,:)=0;
    D2_(:,1:3:szx)=nan;
    D2_(:,3:3:szx)=nan;
    D2_(1:3:szy,:)=nan;
    D2_(3:3:szy,:)=nan;
    
    A = A_;
    A(~isnan(D2_))=0;
    clear imgCol
    O(:,:,1)=D; % red
    O(:,:,2)=A; % green
    O(:,:,3)=Aorig; % blue
    imwrite(O,overlayFN,'WriteMode','append','Compression', 'none'); % 
    waitbar(i/Frames,hw,'processing frames...')
end
close(hw);

mxD = max(D(:));
mxA = max(A(:));








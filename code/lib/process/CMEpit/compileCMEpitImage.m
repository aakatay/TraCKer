% reads pit data from mat files sums and overlaps with timelapse mean
close all;
clear all;

cellType = 'DYN';
%cellType = 'AP2';

pitFN = rdir('**\pitRec.mat'); % recruitment images
tlFN = rdir('**\selPitMean*.tif'); % timelapse mean
% output filename
pitOutFN = 'pitALLRecSum.tif';
overlayOutFN = 'pitRecTLoverlay.tif';

pitREC = [];
nCells = numel(pitFN);
for i = 1:nCells
    fn = pitFN(i).name;
    load(fn);
    pitREC = cat(3,pitREC,pitRec);
end
pitRECsum = uint16(sum(pitREC,3));
npits = size(pitREC,3);
imwrite(pitRECsum,pitOutFN)

% time lapse image
TL = imread(tlFN(1).name);

%size match 
sz = size(pitRECsum,1);
sz2 = size(TL,1);
szpad = (sz2-sz)/2;
pitRECsumPAD = padarray(pitRECsum,[szpad szpad]);

% intensity scale
mxTL = max(TL(:));
mxREC = max(pitRECsum(:));
sc = uint16(round(mxTL/mxREC*2)); 

% overlay
recTLoverlay = uint16(pitRECsumPAD*sc);
recTLoverlay(:,:,2) = TL;
recTLoverlay(:,:,3) = 0;

% filename
overlayOutFN = [overlayOutFN(1:end-4) '-' cellType '_' num2str(npits) 'pits_' num2str(nCells)  'cells.tif'];
imwrite(recTLoverlay,overlayOutFN,'Compression', 'none') 

imagesc(recTLoverlay)
maximize; axis image;
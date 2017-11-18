% using pit crops from (1)time lapse and (2) recruitment data
% finds the average intensity and compare detected edge
clear all

recFN = 'E:\MATLAB\SUM-AP2\170528-SUM\cell7';
pitImgFN = 'E:\MATLAB\SUM-AP2\170521_SUM-AP2\pitselection';
PWD = pwd; 


% pit image
cd(pitImgFN);
pitImg = imread('selPitMean.tif');
mx = double(max(pitImg(:)));

%
cd(recFN);
ixPit = [4 9]; % 14 13 
load pitDetection; % pit rec images
R = imread('binImgRcrtSum.tif');
ws = 9; %windowsize
m =1; % magnification
for i=1:numel(ixPit) % each pit
    
    r = pr(:,:,ixPit(i));
    mxr = max(r(:));
    [Xc,Yc] = centOfMassLocCore(r,16);
    pitRec = pitImg;
    pitRec(9:24,9:24,2) = r/mxr*mx;
    pitRec(:,:,3)=0;
    imwrite(uint16(pitRec),sprintf('pitOverlay%02i_ix%02i.tif',i,ixPit(i)));
end
clear all;
close all;

% input
isBRT = 1; % ow PRB
isCy3 = 1; % ow. GFP
dataFN = 'C:\MATLAB\TraCKer\code\TIRFcalibration\fluorescenceSpectra\spectrumData';
% read input data
[fc,fq,fl,ls,sys] = readSpectraData(isBRT,isCy3,dataFN);

%% display param
disp.FaceAlphaVal = 0.2;
disp.EdgeAlphaVal = 0.7;

% function call
[mnmx] = fluorescenceSpectra(fc,fq,fl,ls,sys,disp,isBRT,isCy3,dataFN); % find the boundaries
fluorescenceSpectra(fc,fq,fl,ls,sys,disp,isBRT,isCy3,dataFN,mnmx); % calculate
    
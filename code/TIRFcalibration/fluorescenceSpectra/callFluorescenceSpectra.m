clear all;
close all;
%% display param
disp.FaceAlphaVal = 0.2;
disp.EdgeAlphaVal = 0.7;

sysBRT.fluoCollect = 0.3; % NA=
sysPRB.fluoCollect = 0.35; % NA=

fcPRB = 'spectrumData\fcPRB.mat';
flCY = 'spectrumData\flCY.mat';
fqBRT = 'spectrumData\fqBRT.mat';
fcBRT = 'spectrumData\fcBRT.mat';
%% load Data
load(fcPRB);    %(fc) filter cube 
load(fcBRT);    %(fc) filter cube 
load(flCY);     %(fl) fluorophores
load(fqBRT);    %(fq) quadview/CSU
lsPRB = [488 561 640];
lsBRT = [473 532 635];
ls.X = lsBRT;   %(ls) lasers
ls.P = [1 1 1]; % laser power


sys = sysBRT;

%% shape input data to structs
% fc
fc.X{1} = fcEx(:,1); % spectra (x-axis)
fc.X{2} = fcEm(:,1);
fc.X{3} = fcBs(:,1);
fc.D{1} = fcEx(:,2); % values (y-axis)
fc.D{2} = fcEm(:,2);
fc.D{3} = fcBs(:,2); 
% fqm
fqm.X{1} = fqm1(:,1); % spectra (x-axis)
fqm.X{2} = fqm2(:,1);
fqm.X{3} = fqm3(:,1);
fqm.D{1} = fqm1(:,2); % values (y-axis)
fqm.D{2} = fqm2(:,2);
fqm.D{3} = fqm3(:,2); 
% fqd
fqd.X{1} = fqd1(:,1); % spectra (x-axis)
fqd.X{2} = fqd2(:,1);
fqd.X{3} = fqd3(:,1);
fqd.D{1} = fqd1(:,2); % values (y-axis)
fqd.D{2} = fqd2(:,2);
fqd.D{3} = fqd3(:,2); 

fq.fqm = fqm; fq.fqd = fqd;

[mnmx] = dispFluorescenceSpectra(fc,fq,fl,ls,disp);

dispFluorescenceSpectra(fc,fq,fl,ls,disp,mnmx);

dispChanIntensity
clear all;
close all;
isBRT = 1; % ow PRB
isCy3 = 1; % ow. GFP
%% display param
disp.FaceAlphaVal = 0.2;
disp.EdgeAlphaVal = 0.7;

sysBRT.fluoCollect = 0.3; % NA=
sysPRB.fluoCollect = 0.35; % NA=

fcBRT = 'spectrumData\fcBRT.mat';
fcPRB = 'spectrumData\fcPRB.mat';
flCY3 = 'spectrumData\flCY.mat';
flGFP = 'spectrumData\flGFP.mat';
fqBRT = 'spectrumData\fqBRT.mat';
fqPRB = 'spectrumData\fqPRB.mat';
lsPRB = [488 561 640];
lsBRT = [473 532 635];
%% load Data
if isBRT
    load(fcBRT);    %(fc) filter cube 
    if isCy3
        load(flCY3);     %(fl) fluorophores
    else
        load(flGFP);     %(fl) fluorophores
    end
    load(fqBRT);    %(fq) quadview/CSU
    ls.X = lsBRT;   %(ls) lasers
    ls.P = [1 1 1]; % laser power
else
    load(fcPRB);    %(fc) filter cube
    if isCy3
        load(flCY3);     %(fl) fluorophores
    else
        load(flGFP);     %(fl) fluorophores
    end
    load(fqPRB);    %(fq) quadview/CSU
    ls.X = lsPRB;   %(ls) lasers 
    ls.P = [1 1 1]; % laser power
end
sys = sysBRT;

%% shape input data to structs
% fc
fc.X{1} = fcEx(:,1); % spectra (x-axis)
fc.X{2} = fcEm(:,1);
fc.X{3} = fcBs(:,1);
fc.D{1} = fcEx(:,2); % values (y-axis)
fc.D{2} = fcEm(:,2);
fc.D{3} = fcBs(:,2); 

if isBRT % quadview
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
else % CSU
    % fqm
    fqm.X{1} = fqm1(:,1); % spectra (x-axis)
    fqm.D{1} = fqm1(:,2); % values (y-axis) 
    fq.fqm = fqm; fq.fqd = [];
    
end

[mnmx] = dispFluorescenceSpectra(fc,fq,fl,ls,disp,isBRT,isCy3);

dispFluorescenceSpectra(fc,fq,fl,ls,disp,isBRT,isCy3,mnmx);

dispChanIntensity
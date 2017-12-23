% generates MAT data files 
% run in spectrumData\

%% KURAL LAB (PRB)
%filter cube (Chroma 89401 - ET - DAPI/FITC/TRITC/CY5 Quad)
fnX = '89101x.txt';
fnM = '89101m.txt';
fnB = '89100bs.txt';
fcEx = readSpectra(fnX);
fcEm = readSpectra(fnM);
fcBs = readSpectra(fnB); 
fcEx = reshape(fcEx,2,numel(fcEx)/2)';
fcEm = reshape(fcEm,2,numel(fcEm)/2)';
fcBs = reshape(fcBs,2,numel(fcBs)/2)';
save('fcPRB','fcEx','fcEm','fcBs');
% CSU
fqm1 = [fcBs(:,1) fcBs(:,2)*0+1]; % none
save('fqPRB','fqm1');

%% FISHEL LAB (BRT)
fnB = 'semrockDi03-R405-488-532-635-t1-25x36.txt';
fcBs = readSpectra(fnB);
fcBs = reshape(fcBs,2,numel(fcBs)/2)';
fcEx = [fcBs(:,1) fcBs(:,2)*0+1]; % none
fcEm = [fcBs(:,1) fcBs(:,2)*0+1]; % none
save('fcBRT','fcEx','fcEm','fcBs');

% quadview
fqd1 = readSpectra('QVd1_538semrock.txt');
fqd1 = reshape(fqd1,2,numel(fqd1)/2)';
fqd2 = readSpectra('QVd2_640semrock.txt');
fqd2 = reshape(fqd2,2,numel(fqd2)/2)';
fqd3 = readSpectra('QVd3_740chroma.txt');
fqd3 = reshape(fqd3,2,numel(fqd3)/2)';

fqm1 = readSpectra('QVm1_511-20-25semrock.txt');
fqm1 = reshape(fqm1,2,numel(fqm1)/2)';
fqm2 = readSpectra('QVm2_593-46semrock.txt');
fqm2 = reshape(fqm2,2,numel(fqm2)/2)';
fqm3 = readSpectra('QVm3_690-50chroma.txt');
fqm3 = reshape(fqm3,2,numel(fqm3)/2)';
save('fqBRT','fqd1','fqd2','fqd3','fqm1','fqm2','fqm3');




%% fluo CY
cy2fn = 'Alexa Fluor 488.csv';
cy3fn = 'Cy3.csv';
cy5fn = 'Cy5.csv';

cy2 = readSpectra2(cy2fn);
cy3 = readSpectra2(cy3fn);
cy5 = readSpectra2(cy5fn);


% spectra (x-axis)
fl.X{1} = cy2(:,1);
fl.X{2} = cy3(:,1);
fl.X{3} = cy5(:,1);

% values (y-axis)
% excitation
fl.Dx{1} = cy2(:,2)/100;
fl.Dx{2} = cy3(:,2)/100;
fl.Dx{3} = cy5(:,2)/100;
% emission
fl.Dm{1} = cy2(:,3)/100;
fl.Dm{2} = cy3(:,3)/100;
fl.Dm{3} = cy5(:,3)/100;
cy3Q = 0.31; % 
cy3e = 150000; % extinction coeff
fl.B = [1 1 1]*1e-6; % fluo brightness
save('flCY','fl');


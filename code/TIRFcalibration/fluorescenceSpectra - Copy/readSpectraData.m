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

%% FISHEL LAB (BRT)
fnB = 'semrockDi03-R405-488-532-635-t1-25x36.txt';
fcBs = readSpectra(fnB);
fcBs = reshape(fcBs,2,numel(fcBs)/2)';
fcEx = nan;
fcEm = nan;
save('fcBRT','fcEx','fcEm','fcBs');

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
save('flCY','fl');



%% display param
disp.FaceAlphaVal = 0.2;
disp.EdgeAlphaVal = 0.7;

fcPRB = 'spectrumData\fcPRB.mat';
flCY = 'spectrumData\flCY.mat';

%% load Data
load(fcPRB); %(fc) filter cube 
load(flCY); % (fl) fluorophores
ls.X = [405 488 561 640]; % (ls) lasers

% shape data
% spectra (x-axis)
fc.X{1} = fcEx(:,1);
fc.X{2} = fcEm(:,1);
fc.X{3} = fcBs(:,1);
% values (y-axis)
fc.D{1} = fcEx(:,2);
fc.D{2} = fcEm(:,2);
fc.D{3} = fcBs(:,2); 

[mnmx] = dispFluorescenceSpectra(fc,fl,ls,disp);

dispFluorescenceSpectra(fc,fl,ls,disp,mnmx);

%% 
    tx = 0:3:300-3;         % Time vector for original signal
    x = sin(2*pi*tx/300);   % Define a sinusoid 
    ty = 0:2:300-2;         % Time vector for resampled signal        
    y = resample(x,ty);    % Change sampling rate
    plot(tx,x,'+-',ty,y,'o:')
    legend('original','resampled');
    xlabel('Time')
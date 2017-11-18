px = 106; % [nm]
sg_ = 77; % [nm] => abbe res. : 183nm
sg = sg_/px; % [px]
%% parameters
gain = 300;
% SNR (aim single photon sensitivity)
noiseRat = 0.5;
SNR = 1/noiseRat;
% emission

excPow = 60; % [mW]
acqTime = 1.9; % [ms]
illumAreaPx = 50*50;
illumArea = illumAreaPx*px^2*1e-6; % um^2
excPowPerArea = excPow/illumArea; % [mW/um^2] = ExcDens (excitation density)
exc_ms = excPowPerArea % excitation in msec.
phPerFluoPerExc_ms = 
phPerSecPerFluo = excPow; % photons per fluo. per frame
phPerFramePerFluo = 20; % photons per fluo. per frame

imgSz = 100;
a = zeros(imgSz,imgSz); % image
x_ = 1:imgSz;
y_ = x_;
[xg,yg]=meshgrid(x_,y_); % x-y grid

Y = [50  20  10+dy  38  34];
X = [50  50  10+dx  34  34];
P = [1 1 1 1 1]*phPerFramePerFluo; % intensity (photons)

for i = 1:
    y=Y(i);x=X(i);p=P(i);
    a = a + p*exp(-((x-xg).^2+(y-yg).^2)/2/sg^2);
end
    
% add background noise
gausBlur = fspecial('gaus',3,2);
noiseBck = rand(size(a));
noiseBck = imfilter(noiseBck,gausBlur,'symmetric');

% add readout noise
noiseRdOut = rand(size(a));


noise = noise*mnP/SNR/mxNoise;

a = noise + a ;

%% fluo collection
NA = 1.49;
RI = 1.52;
T = asin(NA/RI);
collectionRatio = (1-cos(T))/2;

%% surfIllumEffect
surfIllumEffect = 10;
inclineAngle = 


%% 
wl = 0.5e-6; % [m]
cameraDetector = 2^26; % [um^2] 16* 16*512*512
mag = 100*1.5;
FOV = cameraDetector/mag; % [um^2]
LS0 = 1; % [mW] laser source power
filterCube = 1; % filtering after dichroic
cropFOV = 1; % ratio of power in FOV
LS = LS0*filterCube*cropFOV;
LF = LS/FOV*1e8; % laser excitation flux (nW/um2 = mW/mm2)

h = 6.626070040e-34; % [J.s] Planck constant
c = 3e8; % [m.s-2]speed of light
Na = 6.0221409e23; % Avagadro's number

pE = h*c/wl; % [w] 
pF = LF*1e-9/pE; % photon flux (number of photons/um2)

ex = 150000; % [L.mol-1.cm-1]  molar extinction coefficient 
QY = 0.31; % quantum yield
ax = ex/Na; % [cm2/molec] absorption cross section
axNM2= ax*1e14; % [nm2/molec] absorption cross section
ax1D = sqrt(axNM2)*1e1; % Armstrong

pEx = pF*ax*1e8; % exciting photons
pEm = QY*pEx;
pDet0 = pEm*collectionRatio; % detected emitted photons
pDet = pDet0*surfIllumEffect;
pDet

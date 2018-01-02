
frstzeroX   = ; % first zero positiob of the airy distrib
frstzeroKx  = ;
sqX         = ; % asymmetric squeeze factor
sqKx        = ;
nX          = ; % number of diffraction orders
nKx         = ;
frstzero    = [frstzeroX frstzeroKx];
sq          = [sqX sqKx];
n           = [nX nKx];

[x,kx,I] = genBlizzProfile(frstzero,sq,n);

 

%% spatial distribution
% Gaussian
%% angular distribution
% of different illuminations
clear
close all;

szXY    = 100;
sgmXY   = 20;
szZ     = 200;
sgmZ    = 20;

%% 3D Gauss function
x = -szXY:szXY; % grid vectors
z = -szZ:szZ;
[a c]= meshgrid(x,z); % GRID
arg   = -(a.*a)/(2*sgmXY*sgmXY)-c.*c/(2*sgmZ*sgmZ);
G     = exp(arg); % exponential
G  = G/sum(G(:)); % normalize
contour(G); % display Gauss function

%% regenerate 2D slices of the 3D Gauss function
zVal = G(:,szXY+1);
a = x;
for i = 1:length(zVal)
    arg   = -a.*a/(2*sgmXY*sgmXY);
    g(i,:)     = zVal(i)*exp(arg); % exponential
end

figure; imagesc(G-g); colorbar;



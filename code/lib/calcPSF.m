clear;
close all;
goRoot;
CD = cd;
%% get the folder names of DS(raw) and 3Drotated data
if ismac 
    fPathDS = rdir('data/DSdata/**/*DS','isdir==1');
    fPath3D = rdir('data/Rotated3D/**/*dorsal3D','isdir==1');
else
    fPathDS = rdir('data\DSdata\**\*DS','isdir==1');
    fPath3D = rdir('data\Rotated3D\**\*dorsal3D','isdir==1');
end
nData = size(fPathDS,1); % number of datasets
%fRange1 =1; fRange2=nData;
fRange =3; fRange2=3;
for i = fRange1 : fRange2
    %% DS folder (rotAngle & phiTheta)
    goRoot;
    cd(fPathDS(i).name); 
    load lateralTilt; % phiTheta
    tiltAngle(:,i)=phiTheta;
    % read rotation angle
    fAngle = fopen(char('rotAngle.txt'),'r'); 
    angle = fread(fAngle);
    rotAngle(i) = str2num(char(angle'));
    fclose(fAngle);
    cd(CD); % return root
    %% 3D folder 
    cd(fPath3D(i).name); 
    d = rdir('traceData-*');
    [Y,I] = sort([d.datenum]); 
    traceDataPath(i) = {d(I(end)).name}; % load the most recent traceData
    
    display(sprintf('---> %s \n', fPathDS(i).name));
end
%% get trace volume for selected spots
ixTraces(:,1) = 1:10
for i = fRange1 : fRange2
    traceVol = dispTraceVol(traceDataPath(i),rotAngle(i),tiltAngle(:,i))
    
end
% run in a folder with IN and OUT folders
% IN:DS and raw stacks
% reads and measures DS image
clear
close all

% pixel sizes
pxSz = 0.104; % camera's pixel size
stepSz = 0.38; % stage's step distance

imgFileDS = 'IN\DSimageZYX.tif';
imgInfoDS = imfinfo(imgFileDS);
imgFileRAW = 'IN\RAWimageZYX.tif';
imgInfoRAW = imfinfo(imgFileRAW);
nzDS = imgInfoDS.Width;
nyDS = imgInfoDS.Height;
nxDS = numel(imgInfoDS); % layers
nzRAW = imgInfoRAW.Width;
nyRAW = imgInfoRAW.Height;
nxRAW = numel(imgInfoRAW); % layers
disp(sprintf('RAW IMAGE SIZE: \nx:%i,\ny:%i,\nz:%i',nxRAW,nyRAW,nzRAW));
disp(sprintf('DS IMAGE SIZE: \nx:%i,\ny:%i,\nz:%i',nxDS,nyDS,nzDS));
for i = 1: nxDS
    imgDS(:,:,i) = double(imread(imgFileDS,i)); % ZYX
end
%XZY
stackWrite(flipud(permute(imgDS,[2 3 1])),'OUT\DSimageXZY.tif'); % imgDS_XZY
% NOTE : DSimageXZY image seem to be problematic may be DS data is not good

for i = 1: nxRAW
    imgRAW(:,:,i) = double(imread(imgFileRAW,i)); % ZYX
end
%XZY
stackWrite(flipud(permute(imgRAW,[2 3 1])),'OUT\RAWimageXZY.tif'); % imgDS_XZY

%% find the angle of the DS XZ projection
% y projection
for i = 1:nzDS
    for j = 1:nxDS
        imgXZds(j,i) = mean(imgDS(imgDS(:,i,j)>100,i,j));
    end
end
imgXZds = flipud(imgXZds);
imagesc(imgXZds);
imgXZds=imgXZds*2^15/max(imgXZds(:));
imwrite(uint16(imgXZds),'OUT\imgDS_XZproj.tif');


% same for raw data
for i = 1:nzRAW
    for j = 1:nxRAW
        imgXZraw(i,j) = mean(imgRAW(imgRAW(:,i,j)>100,i,j));
    end
end
imgXZraw = flipud(imgXZraw);
imagesc(imgXZraw);
imgXZraw=imgXZraw*2^15/max(imgXZds(:));
imwrite(uint16(imgXZraw),'OUT\imgRAW_XZproj.tif');


% angle is defined as the angle betw.:
% the ZY plane and stage's translation axis (X)
% angle found as 18.5 degrees (-108.5)(see XZ angle_18.5degree.PNG) using nonisotrpic image
% convert to the real angle if found around 50 degrees
thetaReal = atan(tan(18.5/180*pi)*stepSz/pxSz)/pi*180;
thetaReal = atan(tan(15/180*pi)*stepSz/pxSz)/pi*180;
disp(sprintf('real theta found as : %.02f',thetaReal));
% this angle can be misleading as the plate of the cell might not be
% parallel to the x axis translation

% y projection
for i = 1:nxDS
    for j = 1:nzDS
        imgXZds(i,j) = mean(imgDS(imgDS(:,i,j)>0,i,j));
    end
end
imagesc(imgXZds)


%% actual theta angle is found as 55 from shear corrected raw data 
% see XZ angle_isotropic_55degree.PNG


dx = stepSz;
dy = pxSz;
dz = pxSz*sin(35/180*pi);

disp(sprintf('PIXEL SIZES: \ndx:%.02f,dy:%.02f,dz:%.02f,',dx,dy,dz))


r = stepSz/pxSz;
Himg = 262*pxSz;
Hreal = Himg*sqrt(2*r^2/(r^2+1))


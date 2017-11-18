% rotate in Z-axis by 45 deg. ( XYZ -> X'Y'Z)
% save stack in X' Z Y'


clear; close all;

%browse main folder
t = 1;
while 7~=exist('DSrotate') % not in the main folder
    cd('..');
    t = t+1;
    if t>10
        display('change the current directory to DSrotate main folder');
    end
end 
    
CD = cd; % current directory
addpath(strcat(cd,'/DSrotate'));

fMainPathImg3D = 'rotated3D/'; % the folder where rotated DS images are stored
fPath3D = dir(fMainPathImg3D); % folders for each cell
%fPath3D = fPath3D(1:6);
%fPath3D(6).name = '2013_10_30 (Kirchausen) AP2 - leading edge sample_1-4'

img3D=0;
for i = 5 : size(fPath3D,1)
    clear img3D
    cd(char(strcat(fMainPathImg3D,fPath3D(i).name)));
    display(sprintf('---> %s \n', fPath3D(i).name));
    %% read file names
    imgNames = rdir('*.tif');
    imgList = struct2cell(imgNames);
    K=1;
    tifName = imgList(1,1:K);
    % size of the image
    nImage = size(tifName,2);
    imageInfo=imfinfo(char(strcat(tifName(1))));
    numFrames=length(imageInfo);
    imSize=[imageInfo(1).Height,imageInfo(1).Width,numFrames];
    rotAngle = atan(imSize(1)/imSize(2))/pi*180
    % maxProj loop
    %for K = 1:1
        display(sprintf('image : %d/%d \n',K,nImage));
        %for z=1:numFrames % # of frames
         %   img3D(:,:,z) = imread(char(strcat(tifName(K))),z); 
            %img3D = flipdim(img3D,2);
        %end
        if 7 ~= exist(strcat(char(CD),'/RotatedXY'))
            mkdir(char(CD),'/RotatedXY');
        end
        if 7 ~= exist(strcat(CD,'/RotatedXY/',fPath3D(i).name))
            mkdir(strcat(CD,'/RotatedXY/',fPath3D(i).name));
        end
        dFolder = strcat(CD,'/RotatedXY/',fPath3D(i).name);
        fnameRtdXY = strcat(CD,'/maxProj/',fPath3D(i).name,'/',tifName(K));
        rotateXYCore(imgNames(1).name,fPath3D(i).name,CD,dFolder,imSize,rotAngle);
    %end
    cd(CD);    
end
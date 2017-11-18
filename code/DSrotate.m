

ftif = dir('DS\*.tif');


%Extract the number of frames from the tiff file's information. The
%number of columns in structure is equal to the number of frames. 
imageInfo=imfinfo(ftif(1).name);
numFrames=length(imageInfo);
%Size the movie matrix
imSize=[imageInfo(1).Height,imageInfo(1).Width,numFrames];
%nImage = 3;
for K = 1:nImage
    display(sprintf('image : %d/%d \n',K,nImage));
    % read rotation angle
    fAngle = fopen(char(strcat(fPath,'rotAngle.txt')),'r'); 
    rotAngle = fread(fAngle);
    rotAngle = str2num(char(rotAngle'));
    fclose(fAngle);
    %zRange clip
    zRange12=0; % no cropping
    if 2==exist('zRange12.mat') 
        load('zRange12.mat');
    end
    DSrotateCore(tifNameRd(K),tifNameWr(K),fPath,CD,dFolder(i+jj).name,imSize,rotAngle,zRange12,isReConfig);
end

%% call skewCorrInterp
skewCorrInterp

%% Max proj time stacks 
f1_ = strcat(CD,'/data/Rotated3D/',dFolder(i+jj).name,'/3D/maxProj3D.tif');
f1 = rdir(strcat(CD,'/data/Rotated3D/',dFolder(i+jj).name,'/3D'),'strfind(name, ''maxProj-stack'')');
f2_ = strcat(CD,'/data/Rotated3D/',dFolder(i+jj).name,'/dorsal3D/maxProjDorsal3D.tif');
f2 = rdir(strcat(CD,'/data/Rotated3D/',dFolder(i+jj).name,'/dorsal3D'),'strfind(name, ''maxProj-stack'')');
f3_ = strcat(CD,'/data/Rotated3D/',dFolder(i+jj).name,'/ventral3D/maxProjVentral3D.tif');
f3 = rdir(strcat(CD,'/data/Rotated3D/',dFolder(i+jj).name,'/ventral3D'),'strfind(name, ''maxProj-stack'')');

for i = 1 : length(f1)
    if i == 1
        imwrite(imread(char(f1(i).name)),f1_,'tiff');
    else
        imwrite(imread(char(f1(i).name)),f1_,'tiff','WriteMode','append');
    end
    
end
for i = 1 : length(f2)
    if i == 1
        imwrite(imread(char(f2(i).name)),f2_,'tiff');
    end
    imwrite(imread(char(f2(i).name)),f2_,'tiff','WriteMode','append');
end
for i = 1 : length(f3)
    if i == 1
        imwrite(imread(char(f3(i).name)),f3_,'tiff');
    end
    imwrite(imread(char(f3(i).name)),f3_,'tiff','WriteMode','append');
end
return;
% delete the temp max proj files
f = rdir('maxProj-stack*');
delete(f.name);
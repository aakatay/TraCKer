

ftif = dir('DS\*.tif');
%Extract the number of frames from the tiff file's information. The
%number of columns in structure is equal to the number of frames. 
nImage=length(ftif);
imageInfo=imfinfo(['DS\' ftif(1).name]);
numFrames=length(imageInfo);
%Size the movie matrix
imgSize=[imageInfo(1).Height,imageInfo(1).Width,numFrames];
%nImage = 3;


% find the shearing angle :thISO 
% th : shearing of the acquired image (noniso sampled)
pxSzX = 0.104; %pixel sizes
pxSzZ = 0.38;
th = 90-32.8;
pxZvsX = pxSzZ / pxSzX; 

isReConfig = 0;
nImage = 2;
for K = 1:nImage
    display(sprintf('image : %d/%d ',K,nImage));
    % read rotation angle
    fAngle = fopen('rotAngle.txt','r'); 
    rotAngle = fread(fAngle);
    rotAngle = str2num(char(rotAngle'));
    fclose(fAngle);
    %zRange clip
    tifNameRd = strcat('DS\',ftif(K).name);
    tifNameWr = sprintf('stack_%03i.tif',K);
    DSrotateCore(tifNameRd,tifNameWr,imgSize,rotAngle,pxZvsX,isReConfig);
end

%% call skewCorrInterp
%skewCorrInterp

%% Max proj time stacks 
f1_ = 'Rotated3D\3D\maxProj3D.tif';
f1 = rdir('Rotated3D\3D','strfind(name, ''maxProj-stack'')');
f2_ = 'Rotated3D\dorsal3D\maxProjDorsal3D.tif';
f2 = rdir('Rotated3D\dorsal3D','strfind(name, ''maxProj-stack'')');
f3_ = 'Rotated3D\ventral3D\maxProjVentral3D.tif';
f3 = rdir('Rotated3D\ventral3D','strfind(name, ''maxProj-stack'')');

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
    else
        imwrite(imread(char(f2(i).name)),f2_,'tiff','WriteMode','append');
    end
end
for i = 1 : length(f3)
    if i == 1
        imwrite(imread(char(f3(i).name)),f3_,'tiff');
    else
        imwrite(imread(char(f3(i).name)),f3_,'tiff','WriteMode','append');
    end
end
% delete the temp max proj files
delete(f1.name)
delete(f2.name)
delete(f3.name)
return;
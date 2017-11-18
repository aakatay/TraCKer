% extracts the high layers from the data to determine the outline of the
% cell
%clear;
nZ = 25; % # of z planes
fname = 'C2-july18_acini from july 10 imaged after 28  hours_Subset 1 red channel.tif';
imageInfo=imfinfo(fname);
numFrames=length(imageInfo);
%Size the movie matrix
imSize=[imageInfo(1).Height,imageInfo(1).Width,numFrames];

filenameMaxProj_out = 'maxProj.tif';
filenameStack_out = 'stack_';
filenameMaxProj_out = 'maxProj.tif';

for i = 1 : numFrames
    A(:,:,i) = imread(fname,i);
end
isCrop = 0;
yy = 1:250;
xx = 1:250;
h = waitbar(0,'slicing data stack...');
for i = 1 : size(A,3)/nZ
    for j = 1 : nZ
        if isCrop, temp = A(yy,xx,(i-1)*nZ+j); else, temp = A(1:end,1:end,(i-1)*nZ+j); end;
        %Zstack(:,:,j) = temp(200:400,150:350);
        Zstack(:,:,j) = temp;
    end
    fn = sprintf('%s%02i.tif', filenameStack_out,i);
    stackWrite(Zstack,fn); 
    maxProj(:,:,i) = max(Zstack,[],3);
    waitbar(i / size(A,1)*nZ)
end
stackWrite(maxProj,filenameMaxProj_out);
close(h)

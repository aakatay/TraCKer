% extracts the high layers from the data to determine the outline of the
% cell
%clear;
close all;
fname = 'C1-june29 acini 10 days old image 1.tif';
imageInfo=imfinfo(fname);
numFrames=length(imageInfo);
%Size the movie matrix
imSize=[imageInfo(1).Height,imageInfo(1).Width,numFrames];

filenameMaxProj_out = 'maxProj.tif';
nZ = 36;
if ~exist('A')
    for i = 1 : numFrames
        A(:,:,i) = imread(fname,i);
    end
end

Zstack = ones(imSize(1),imSize(2));
maxProj = Zstack;
h = waitbar(0,'slicing data stack...');
for i = 1 : size(A,3)/nZ
    for j = 1 : nZ
        Zstack(:,:,j) = A(:,:,(i-1)*nZ+j);
    end
    maxProj(:,:,i) = max(Zstack,[],3);
    waitbar(i / size(A,1)*nZ)
end
stackWrite(maxProj,filenameMaxProj_out);
close(h)

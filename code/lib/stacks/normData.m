if exist('isCallRunBaLM')
    clear all;
end

fname = 'data.tif';
fOut = 'dataNorm.tif';

imageInfo=imfinfo(fname);
numFrames=length(imageInfo);
%Size the movie matrix
imSize=[imageInfo(1).Height,imageInfo(1).Width,numFrames];
for i = 1 : numFrames
    A(:,:,i) = double(imread(fname,i));
end

x1 = 1;
y1 = 1;
A = A(x1:end,y1:end,:);

A = A-min(A(:));
maxData = max(A(:));

A = A/maxData*(2^16-1);

stackWrite(uint16(A),fOut);
clear
if isFile
    fname = 'NBINs.tif';
    imageInfo=imfinfo(fname);
    numFrames=length(imageInfo);
    %Size the movie matrix
    imSize=[imageInfo(1).Height,imageInfo(1).Width,numFrames];
    xx = imSize(1);
    yy = imSize(2);
    for i = 1 : numFrames
        A(:,:,i) = double(imread(fname,i));
    end
else 
    matFile = dir('nbins*.mat');
    load(matFile.name);
end
% select frames
f1 = 1; f2 = 300;

% conv. 
convRect = ones(1,10);
for i = 1:xx
    for j = 1:yy
        Aconv(i,j,:)=conv(A(i,j,f1:f2),convRect,'same');
        Arecruit = imregionalmax(Aconv);
    end
end

for i = 1:size(Aconv,3)
    img = Aconv(:,:,i);
    
end

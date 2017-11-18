fname = 'maxProjDorsal3D.tif'
imageInfo = imfinfo(fname);
Frames = length(imageInfo);

convFunc = ones(7,5);
convFunc = convFunc/numel(convFunc);

for i = 1:Frames
    img = double(imread(fname,i));
    imgFilt = img - conv2(img,convFunc,'same')
    imagesc(imgFilt)
    imagesc(img)
end
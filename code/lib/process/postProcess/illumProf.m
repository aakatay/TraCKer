% defines a valid intensity level for bleaching to be used as a mask for
% selecting structures
clear all;
close all;

illumFNdir = rdir('..\..\','\*illum_frm*');
illumFNdir = rdir('..\..\**\*','strfind(name,''illum_'')');
cropXY = 

illumMeanFN = 'illumAll.tif';
m = numel(illumFNdir)
k=1;
for j = 1:m
    illumFN=illumFNdir(j).name
    iminf = imfinfo(illumFNdir(j).name);
    n=numel(iminf);
    for i = 1:n
        A(:,:,k)=imread(illumFN,i);
        k=k+1;
    end
end

I = mean(A,3);
imagesc(I); axis image
imwrite(uint16(I),illumMeanFN);
Icv = conv2(I,ones(10),'same');
figure(2)
imagesc(Icv); axis image
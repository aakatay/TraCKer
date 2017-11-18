
clear all;close all
load Spots % load spot images
nSpot = size(S,4);
nAcq = 19;
nFrames = 95;

% select spots
sp1=1; sp2=nSpot;
nAcq1 = 1;nAcq2 = nAcq;
nAcq1 = 1; nAcq2 = 2;
frm1= 1; frm2= nFrames;
frm1= 1; frm2= 2;

% plot
for i = nAcq1:nAcq2
    for j = sp1:sp2
        for k = frm1:frm2
            IMGin = S(:,:,(i-1)*nFrames+k,j);
            [IMG,DATA] = Gaussian_TIFF(IMGin);
            sizex = size(S,1);
            sizey = sizex;
            [X,Y] = meshgrid(double(1:sizex),double(1:sizey));

            plotting(IMGin,DATA)
            waitforbuttonpress
        end
    end
end
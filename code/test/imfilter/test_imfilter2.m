clear;
%% output
fnameFlt = ; % Tiff file 
%% input
load img3D %uint8
img3D = img3D
%% parameters
sz=5; %size
sg=20; %sigma
szLap = 4; % half length of the laplace filter
gaus = gauss1D(sz,sg);
lap=-ones(szLap*2+1,1); lap(szLap+1)=szLap*2;
hOld = size(img3D,3);

%% filter
for i = 1:size(img3D,1)
    for j = 1:size(img3D,2)
        cln = squeeze(img3D(i,j,:)); % column (z-axis)
        cln = conv(double(cln),gaus);
        cln = uint8(conv(cln,lap));
        img3Dflt(i,j,:) = cln;
        [imgShrp(i,j), imgTpg(i,j)] = max(cln); % topography
    end
    h = waitbar(i/size(img3D,1));
end
close(h)

%% size adjustment
hNew = size(img3Dflt,3);
cropSz = (hNew-hOld)/2;
img3Dflt = img3Dflt(:,:,cropSz+1:end-cropSz); % cropping after conv
imgTpg = imgTpg-cropSz; % correct the heigh value after conv

%% display
figure(1)
imagesc(imgTpg);

%% save figures
save(sprintf('imgFlt-sz%d-sg%d-szLap%d.mat',sz,sg,szLap),'imgFlt');
% save stack
imwrite(uint16(img3Dflt(:,:,1)),char(fnameFlt),'tiff');
imwrite(uint16(img3Dflt(:,:,1)),char(fnameRtd),'tiff');
for z = 2:size(img3Dflt,3)
    imwrite(uint16(img3D(:,:,z)),char(strcat(fnameRtd)),'tiff','WriteMode','append');
end
return;
figure(2)
sliceomatic(img3Dflt);


     
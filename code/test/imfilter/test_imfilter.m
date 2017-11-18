load img3D %uint8
sz=20; %size
sg=50; %sigma
gaus = gauss3D(sz,sg);
imgFlt = imfilter(flipdim(img3D,3),gaus,'symmetric');
save(sprintf('imgFlt-sz%d-sg%d.mat',sz,sg),'imgFlt');
imgTpg = ; % topography
sliceomatic(imgFlt);


     
% thresholding
% (IMGfilt,IMG,CoeffThresh) --> (BW,din)
dataFiltDiv = IMGfilt./CoeffThresh;
bin = im2bw(dataFiltDiv,1);
din = uint16(bin).*uint16(IMG);
BW = imregionalmax(din, 8);
if exist('debug.txt')
    ccc=2;
end
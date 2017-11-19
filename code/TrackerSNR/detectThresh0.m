% detects threshold for tracking
% output Coeff.mat

% UPDATE required : requires histogram fit
% E:\MATLAB\TIRFcalibration\data\Ata01_3autoCoeff


% image =======
% filtered image =======
% filtered image peak MASK 
% image PEAKs =======
% histogram peak (image PEAKs)
clear all; 
close all;

dbg = 1;
fname = 'IMG.tif';


% coefficient multiplier
bc = 4; % [nM] background concentration
cm1 = bc/10;
cm2 = bc/27;
coeffMul1 = 1+0.4; 
coeffMul2 = 1-0.15; 
coeffMul1 = 1+cm1; % corrects shift of histogram peak to the right
coeffMul2 = 1-cm2; % corrects for increased intensity fluctuations

coeffMul = coeffMul1*coeffMul2; % 4nM background


imginf = imfinfo(fname);
w = imginf(1).Width;
h = imginf(1).Height;

%% image
IMG = imread(fname,1);
            IMG(90:94,1:4)=IMG(95:99,1:4); % ======remove
[~,c] = loghist(IMG,30,0);
hstep = c(2)-c(1); % hist step

if dbg
    nc = 6; % # columns
    figure(505);
    subplot(2,nc,1)
    imagesc(invert(IMG));
    colormap('gray');axis image; title('image')
    
    subplot(2,6,nc+1)
    [~,c] = loghist(IMG,30,1);
end

%% filtered image 
gausKernelSz = 3;
gausKernelSg = 0.7;
gausKernelSzPadded = 5;

pdSz = (gausKernelSzPadded - gausKernelSz)/2; % pad size
gausHat1=fspecial('gaussian', gausKernelSz, gausKernelSg);
lap=[-1,-1,-1;-1,8,-1;-1,-1,-1];

IMGfilt2 = nan(size(IMG));
IMGfilt1 = imfilter(IMG,gausHat1,'symmetric');
IMGfilt1 = imfilter(IMGfilt1,lap,'symmetric');
IMGfilt = max([IMGfilt1(:) IMGfilt2(:)],[],2);
IMGfilt = reshape(IMGfilt,size(IMGfilt1));
if dbg
    subplot(2,nc,2)
    imagesc(invert(IMGfilt))
    colormap('gray');axis image;title('filtered image')
    imwrite(IMGfilt,'IMGfilt.tif')
    subplot(2,6,nc+2)
end

[~,c] = loghist(IMGfilt,30,dbg);

%% filtered image peak MASK 
BW = imregionalmax(IMGfilt, 8);

%% image PEAKs
IMGpeak = uint16(BW).*IMG;
if dbg
    subplot(2,nc,3)
    imagesc(invert(IMGpeak))
    colormap('gray');axis image; title('image PEAKs')
    imwrite(IMGpeak,'IMGpeak.tif')
    
    subplot(2,6,nc+3)
    [hc,c] = loghist(IMGpeak,100,dbg); % hist count and centers
    xlabel('intensity'); ylabel('counts')
end


%%  histogram peak(image PEAKs)
[pks,locs] = findpeaks(hc);
peakInt1 = c(locs(1));

Coeff = peakInt1*coeffMul;
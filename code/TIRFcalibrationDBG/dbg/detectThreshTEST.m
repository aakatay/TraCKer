% detects threshold for tracking
% output Coeff.mat

% UPDATE required : requires histogram fit
% E:\MATLAB\TIRFcalibration\data\Ata01_3autoCoeff


% image =======
% filtered image =======
% filtered image peak MASK 
% image PEAKs =======
% threshold#1 (image PEAKs)
% image peak MASK  =======
% filtered image PEAKs thresholded with image peak MASK =======
% threshold#2 (filtered image) =======

clear all; 
close all;

dbg = 1;
fname = 'IMG.tif';
coeff1 = 1.5;
coeff1 = [1.5 1.7 1.9 2 2.1 2.2 2.3];
coeff1 = [1.2:0.1:1.7];
coeff1 = [1.4:0.05:2.1];
%coeff1 = 1.5;
coeff2 = 1.5;

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


%% threshold#1 (image PEAKs)
[pks,locs] = findpeaks(hc);
peakInt1 = c(locs(1));
if dbg, subplot(2,6,nc+3);title(sprintf('peak: %i',peakInt1)); end;

legendCoeff1txt = {};
for i = 1:numel(coeff1)
    peakThresh = round(peakInt1*coeff1(i)); 
    legendCoeff1txt = [legendCoeff1txt {  sprintf('%i(x%.01f)',peakThresh,coeff1(i)) }];
    
    %% image peak MASK 
    IMGpeakMASK = IMGpeak;
    IMGpeakMASK(IMGpeak<peakThresh) = 0;

    if dbg
        subplot(2,nc,4)
        imagesc(invert(IMGpeakMASK))
        colormap('gray');axis image;title('image peak MASK')
        imwrite(IMGpeakMASK,'IMGpeakMASK.tif')

        subplot(2,6,nc+4)
        if i ==1, 
            [hc,c] = loghist(IMGpeak,100,dbg); % hist count and centers
            xlabel('intensity'); ylabel('counts')
            title(sprintf('peak: %i',peakInt1))
        end
        hold on; line([peakThresh peakThresh],[1 hc(1)]); hold off;
    end

    %% filtered image PEAKs thresholded with image peak MASK
    IMGpeakMASKbw = uint16(im2bw(double(IMGpeakMASK),1));
    IMGfiltMASK = IMGfilt.*IMGpeakMASKbw;

    if dbg
        subplot(2,nc,5)
        imagesc(invert(IMGfiltMASK))
        colormap('gray');axis image;title('filtered image PEAKs')
        imwrite(IMGfiltMASK,'IMGfiltMASK.tif')
        %subplot(2,6,nc+5)
        figure(600)
    end
    if i~=1, hold on; end;
    [hc,c] = loghist(IMGfiltMASK,30,1); % hist count and centers
    hold off;
    figure(505)
    % ================= requires histogram fit =================

    %% threshold#2 (filtered image)
    [pks,locs] = findpeaks(hc);
    peakInt2 = c(locs(1));
    Coeff = round(peakInt2/coeff2); % threshold
    save Coeff Coeff;
    titleTXT = sprintf('coeff:%i',Coeff);

end

subplot(2,6,nc+4)
legend(legendCoeff1txt)
maximize;


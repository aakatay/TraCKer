% tests PSF detection


imgPSFfn=dir('PSF*.tif');

    imgPSF_ = imread(imgPSFfn(1).name);
    k=6;
    imgSpot = imgPSF_(10-k:25+k,20-k:40+k);
    %imgPSF_ = imgPSF_(:,1:40);
    imgPSFbckgrnd = imgPSF_(1:20,1:20);
    bckThrshld = max(imgPSFbckgrnd(:))*1.2; % background threshold
    imgPSF = imgPSF_;
    imgPSF(imgPSF<bckThrshld)=bckThrshld;
    imgPSF = imgPSF-bckThrshld;
    
    %% 1
    subplot(2,3,1)
    imagesc(imgPSF_)
    axis image
    
    %% 2
    subplot(2,3,2)
    imagesc(imgSpot)
    axis image
    title('spot')
    
    %% 3
    subplot(2,3,3)
    imgPSFcrop = imgPSF(10:25,20:40);
    imagesc(imgPSFcrop)
    axis image
    title('PSF')
    
    %% 4 filter
    filtPSF = double(imgPSFcrop/max(imgPSFcrop(:)));
    imgPSFfilt = imfilter(double(imgSpot),filtPSF,'replicate');
    %imgPSFfilt = deconv(double(imgPSF_),filtPSF)
    subplot(2,3,4)
    imagesc(imgPSFfilt)
    axis image
    title('filtered spot')
    
    %% 5 select peaks with a Coeff
    BACKmean=[min(mean(imgSpot,1)),min(mean(imgSpot,2))];
    BACK = min(BACKmean);
    Coeff = max(imgPSFfilt(:)-BACK)/1.2;
    
    %Coeff = 5000;
    imgPSFfiltThresh = imgPSFfilt;
    imgPSFfiltThresh(imgPSFfiltThresh < Coeff) = 0;
    subplot(2,3,5)
    imagesc(imgPSFfiltThresh)
    axis image
    title('thresholded')
    
    %% 6 reginalmax detect
    
    spotDetect = imregionalmax(imgPSFfiltThresh,8);
    subplot(2,3,6)
    spotPeak = (spotDetect).*imgPSFfilt;
    imagesc(spotPeak)
    axis image
    title('local max')
    
    %% 7 localization
    [B,L] = bwboundaries(spotDetect,'noholes');
    
    
    maximize;
    
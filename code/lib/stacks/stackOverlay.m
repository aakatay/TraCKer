function stackOverlay(imgAfn,imgBfn,outFN,offset,mx)
% overlays 2 tiff stacks
% offset : # of frames shift between them
% mx : intensity scale
% imgAfn : RED
% imgBfn : GREEN
% stackOverlay('cy3_4nmBckgrnd_250ms_396pos_tilt8_shft3.5_002.tif','SNRmovieCV.tif','SNRoverlay.tif',[],[14583 20533])
% stackOverlay('cy3_4nmBckgrnd_250ms_396pos_tilt8_shft3.5_001.tif','SNRmovie.tif','SNRoverlay.tif',[],[14583 20533])

    if numel(offset) == 1 & offset == 0, offset = [0 0]; end

    if isempty(offset), offset=zeros(1,2); end;

    if exist(outFN), delete(outFN); end

    infoA = imfinfo(imgAfn);
    infoB = imfinfo(imgBfn);
    nA = numel(infoA)-offset(1);
    nB = numel(infoB)-offset(2);
    n = min(nA,nB);
    
    hw =  waitbar(0);
    for i = 1:n
        imgA = imread(imgAfn,i+offset(1));
        imgB = imread(imgBfn,i+offset(2));
        waitbar(i/nA,hw,'processing frames...')
        imgColor = imgA;
        imgColor(:,:,2) = imgB;
        imgColor(:,:,3) = 0;
        
        imwrite(uint16(imgColor),outFN,'tif','WriteMode','append','Compression', 'none');
    
    end
    close(hw);
end
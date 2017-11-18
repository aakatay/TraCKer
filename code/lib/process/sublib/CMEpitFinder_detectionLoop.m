function [din]   = CMEpitFinder_detectionLoop(IMG,CoeffThresh)
% applies coeff 
dbg =0;

    inFocus_w_1X = 0;
    gausKernelSz = 3;
    gausKernelSg = 0.7;
    gausKernelSzPadded = 5;
    gausHat1=fspecial('gaussian', gausKernelSz, gausKernelSg);
    elev=1/(gausKernelSzPadded^2-gausKernelSz^2);        
    gausHat1 = gausHat1 + elev;
    lap=[-1,-1,-1;-1,8,-1;-1,-1,-1];

    detectFilter; % IMG --> IMGfilt
    detectThreshold; % (IMGfilt,CoeffThresh) --> (BW,din)
    [~,L] = bwboundaries(BW,'noholes');

    if dbg,
        figure(131);
        imagesc(L)
        figure(141)
        imagesc(BW)
        figure(142)
        imagesc(din)
    end
    cc=3;
end
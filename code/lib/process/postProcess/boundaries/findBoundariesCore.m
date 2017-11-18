function [Bt0,Lt0]=findBoundariesCore(R,cvSz)

    %% (1){convolve + fill inside}  (2){invert and convolve} (3){invert back and find boundary}
    % convolve and remove internal space
    Rcv0 = conv2(R,ones(cvSz),'same');
    Rcv = im2bw(Rcv0,0);
    invRcv = 1-Rcv;
    [~,Lb]=bwboundaries(invRcv,8);
    Rcv2 = Lb - 1;
    Rcv2(Rcv2>0)= 1;

    % invert and convolve 
    Rcv3 = im2bw(conv2((1-Rcv2),ones(cvSz),'same'),0);

    % invert back and find boundaries
    [Bt0,Lt0]=bwboundaries(1-Rcv3,'noholes'); % tight boundaries
    
    if 0 
        figure(1000); imagesc(R); axis image;
        figure(1001); imagesc(Rcv0); axis image;
        figure(1002); imagesc(Rcv); axis image;
        figure(1003); imagesc(Lb); axis image;
        figure(1004); imagesc(Rcv2); axis image;
        figure(1006); imagesc(Rcv3); axis image;
    end
end
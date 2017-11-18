
nFrames = 95; % number of frames
nAcq = 19; % number of acquistion sets
nSpot = 8;

load Spots;
for j=1:nSpot
    for i =1:nAcq
        IMGin(:,:,i) = mean(S(:,:,(i-1)*nFrames+[1:nFrames],j),3);
    end
    stackWrite(IMGin,sprintf('spot_%i.tif',j))
    
end

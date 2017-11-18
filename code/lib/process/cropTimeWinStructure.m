% process the recruitment data to generate structures with defined time
% windows
function R = cropTimeWinStructure(binFN,ps,Lt)

    if ~exist(binFN), binFN = ['../' binFN]; end
    nf = numel(imfinfo(binFN));
    for i = 1:nf % each frame
        binImg(:,:,i) = imread(binFN,i);
        binImgTime(:,:,i) = im2bw(binImg(:,:,i),0)*i;
    end
    
    
    ns = size(ps,1);
    R = zeros(size(binImg,1),size(binImg,2));
    for i = 1:ns % each structure
        lt = uint16(Lt(:,:,i));
        f1 = ps(i,4); % time window
        f2 = ps(i,5);
        
        r = sum(repmat(lt,[1 1 f2-f1+1]).*binImg(:,:,f1:f2),3);
        R = R+r;
        
    end

end

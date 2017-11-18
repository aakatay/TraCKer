% joins the output of recruitmentTrack files acc. to their crop coors.
% assumes 2x2 cropping
% assume single overlap size of crops
% directory names
binFN = rdir('*\*binImgRcrtSum.tif');
if numel(binFN)~=4, error('works for 4 crops'), end
for i = 1:4
    binfn = binFN(i).name;
    nmCrop = binfn(1:find(binfn=='\')-1); % directory name
    
    ix_ = find(nmCrop=='_');
    if numel(ix_)==2
        nmCrop = nmCrop(1:ix_(2)-1);
        ix_ = ix_(1);
    end
        
    ixX = find(nmCrop=='X');
    ixY = find(nmCrop=='Y');
    ixx = find(nmCrop=='x');

    xCr(i) = str2num(nmCrop(ix_+1:ixX-1));
    yCr(i) = str2num(nmCrop(ixX+1:ixY-1));
    szXcr(i) = str2num(nmCrop(ixY+1:ixx-1));
    szYcr(i) = str2num(nmCrop(ixx+1:end));
end
    [~,ixXs] = sort(xCr);
    [~,ixYLs] = sort(yCr(ixXs(1:2))); % left
    ixTL = ixXs(ixYLs(1)); % topleft
    ixBL = ixXs(ixYLs(2)); % bottomleft
    [~,ixYRs] = sort(yCr(ixXs(3:4))); % right
    ixTR = ixXs(ixYRs(1)+2); % topright
    ixBR = ixXs(ixYRs(2)+2); % bottomright
    szOvrlp = xCr(ixTL)+szXcr(ixTL) - xCr(ixTR)-1; % overlap size of crops
    m = floor(szOvrlp/2)*4; % margin1
    n = ceil(szOvrlp/2)*4; % margin2
    
    imgTL = imread(binFN(ixTL).name);
    imgTR = imread(binFN(ixTR).name);
    imgBL = imread(binFN(ixBL).name);
    imgBR = imread(binFN(ixBR).name);
    
    imgTL = imgTL(1:end-m,1:end-m);
    imgTR = imgTR(1:end-m,n+1:end);
    imgBL = imgBL(n+1:end,1:end-m);
    imgBR = imgBR(n+1:end,n+1:end);
    
    
    img = [imgTL imgTR; imgBL imgBR];
    dirname = sprintf('join_%iX%iY%ix%i',xCr(ixTL),yCr(ixTL),size(img,1),size(img,2));
    mkdir(dirname);
    imwrite(img,[dirname '\binImgRcrtSum.tif']);

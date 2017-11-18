% using a polygonal mask data crops the data 

if exist('maxProj0.tif')
    fileread = 'maxProj0.tif';
else
    fileread = 'maxProj.tif';
end

load outlineMask
% read image data
for j=1:size(BW,3)
    J(:,:,j)=imread(fileread,j);
end

for i = 1:size(BW,3)
    Jc(:,:,i) = J(:,:,i).*uint16(BW(:,:,i));
end


if strcmp(fileread,'maxProj.tif')
    movefile('maxProj.tif','maxProj0.tif');
end

for t = 1:size(Jc,3) % # of frames  
    if t == 1
        imwrite(uint16(Jc(:,:,t)),'maxProj.tif');
    else
        imwrite(uint16(Jc(:,:,t)),'maxProj.tif','tiff','WriteMode','append');
    end
end
   
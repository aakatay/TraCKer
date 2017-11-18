% test
% see walkAvDiff
%function [stackOut] = walkAv(fname)
    % calculates and saves walking average
    infImg = imfinfo(fname);
    Frames = length(infImg);
    for i = 1:Frames
        img(:,:,i) = imread(fname,i);
        if rem(i,5000) == 0 
            disp(sprintf('frame number %i',i));
        end
    end
    
    WA = 2;
    nEl = size(img,1)*size(img,2);
    img=reshape(img,nEl,size(img,3));
    for i=1:nEl
        imgWA(:,i) = conv(img(i,:),WA,'same');
    end
    imgWA=reshape(imgWA,size(img,1),size(img,2),size(img,3));
    
%end
    


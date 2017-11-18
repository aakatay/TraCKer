function img = correctLateralTilt(img,phiTheta)
%
% correct for basolateral tilt angle
    s1 = size(img);
    
    phi     = phiTheta(1);
    theta   = phiTheta(2);
    %figure(2); imagesc(img(:,:,5))
    
    img = imrotate(img,phi,'bicubic');
    
    img = permute(img,[2 3 1]);
    img = imrotate(img,theta,'bicubic');
    img = permute(img,[3 1 2]);
    %figure(4); imagesc(img(:,:,5))
    img = imrotate(img,-phi,'bicubic');
    %figure(5); imagesc(img(:,:,11))
    s2 = size(img);
    dx = floor((s2(1)-s1(1))/2);
    dy = floor((s2(2)-s1(2))/2);
    img = img(dx:end-dx+1,dy:end-dy+1,:);
    %figure(6); imagesc(img(:,:,11))
    
end


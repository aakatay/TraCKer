function rotAng = calcRotAng(img)
% calculates the rotation angle to be used in rotation of the deskewed
% image
% img : kymograph of the deskewed image


    imgBW = im2bw(img); % threshold (binary image)
    [y,x]=size(img);
    X = imgBW(1,:);
    Y = imgBW(:,end);
    save('rotAng.mat','imgBW','X','Y')

    rotAng = double( -atan( (x-sum(X))/(sum(Y)) )  /pi*180 );
    %rotAng = -17;
end
function rotateXYCore(fName,fPath,CD,dFolder,imgSize,rotAngle)
% rotate the 3D images
    
    %Pre-size the movie matrix
    imgSize = [imgSize(1) imgSize(2) imgSize(3)];
    %imgSize = [2,21,10];

    for z=1:imgSize(3) % # of frames
        img3D(:,:,z) = imread(fName,z); 
        %img3D = flipdim(img3D,2);
    end
    %img3D = permute(img3D,[3 1 2]);
    imgRtd = imrotate(img3D,rotAngle,'bicubic');

    % write to file frame by frame
    fnameRtdXY = strcat(dFolder,'/',fName);
    writeTiffStack(imgRtd,fnameRtdXY,0);
end
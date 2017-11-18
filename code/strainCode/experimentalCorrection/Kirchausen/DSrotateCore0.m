function DSrotateCore(fName,fPath,imgSize,rotAngle,isCalcRotAngle)
% determine the rotation angle    
    A = imread(char(strcat(fPath,'/',fName)));
    %Pre-size the movie matrix
    imgSize = [size(A) 80];
    %imgSize = [2,21,10];

    for x=1:imgSize(1)
        for z=1:imgSize(3) % # of frames
            img2D(:,z,x) = imread(char(strcat(fPath,'/',fName)),z,'PixelRegion',{[x,x],[1,imgSize(2)]});
        end
        if rotAngle == 0 && x <= 5
            % calculate rotation angle
            rotAng = calcRotAng(img2D(:,:,x))
            %rotAng = -17;
            if x == 1
                %pre-define the rotated image
                imgTemp = ones(imgSize(2),imgSize(3));
                imgTempRtd = imrotate(imgTemp,rotAng,'bilinear');
                imgRtdSize = size(imgTempRtd);
                imgRtd = ones(imgRtdSize(1),imgRtdSize(2),imgSize(1)); % rotated image
            end
            % input rotation angle
            if x == 5
                rotAngle = input('enter the rotation angle: ');
                fAngle = fopen(char(strcat(fPath,'/rotAngle.txt')),'w'); 
                fwrite(fAngle,num2str(rotAngle));
                fclose(fAngle)
                return
            end
        end
    end
end
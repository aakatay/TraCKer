clear;
displayHelp = 1;


if displayHelp
    A = 10*ones(10,10);
    A(2:4,2:4) = 22;    % maxima 12 higher than surrounding pixels
    A(6:8,6:8) = 33;    % maxima 23 higher than surrounding pixels
    %A(2,7) = 44;
    a=0;
    A(3,8-a) = 45.2;
    A(4,9-a) = 44.8
    regmax = imregionalmax(A,8)
    subplot(1,2,1)
    imagesc(A)
    subplot(1,2,2)
    imagesc(regmax)
else
    Frames = 26;
    convSz = 2; % extension of the spot size by convolution
    if ~exist('img3D.mat')
        gaus=fspecial('gaussian', 5, 1.5);
        lap=[-1,-1,-1;-1,8,-1;-1,-1,-1];
        reply='stack_001.tif';
        reply='102x108xy_15ms_100P_296G_1.5x_CROP_105.tif';
        for j=1:Frames
            img3D(:,:,j) = imread(reply,j);
            img3D(:,:,j) = imfilter(img3D(:,:,j),gaus,'symmetric');
            img3D(:,:,j) = imfilter(img3D(:,:,j),lap,'symmetric');
        end
        save img3D
    end

    load img3D
    img3D=img3D/100;
    imgSpts = imregionalmax(img3D,18); % spots
    sz = 2*convSz+1;
    convMtrx = ones(sz,sz,sz);
    sliceomatic(double(imgSpts))
    imgSpts = convn(imgSpts,convMtrx);
    imgSpts = logical(imgSpts);
    imgSpts=1-imgSpts;
    img3Dfilt = img3D.*uint16(imgSpts(convSz+1:end-convSz,convSz+1:end-convSz,convSz+1:end-convSz));
    sliceomatic(img3Dfilt);
end


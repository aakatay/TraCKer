    function overlayIMGandDET

        img = imread('regressSpotDetect_001.tif');
        fn_spotWin=dir('spotWin*.mat');
        load(fn_spotWin.name);

        det = zeros(size(img));
        for i=1:size(xBW,1)
            if yBW(i)*xBW(i) == 0 
                continue
            end
            det(yBW(i),xBW(i)) = 5000;
        end
        det = uint16(det);
        img_overlay = img + det;
        subplot(2,2,1)
        imagesc(img); title('image');
        subplot(2,2,2)
        imagesc(det); title('detected');
        subplot(2,2,3)
        imagesc(img_overlay); title('overlay');
        copyfile(fn_spotWin.name,['..\detectedSpots.mat' ]);
        cd('..\')
    end

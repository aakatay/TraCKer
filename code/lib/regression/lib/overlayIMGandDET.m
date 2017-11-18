    function overlayIMGandDET
        try
            img = imread('..\regressSpotDetect_10ms_001.tif');
        catch
            img = imread('regressSpotDetect_10ms_001.tif');
        end
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
        hFig = figure;
        
        set(hFig,'Color',[0 0 0])
        img_overlay = img + det;
        subplot(2,2,1)
        imagesc(img); title('image');
        set(gca,'Parent',hFig,'YColor',[1 1 1]*0.7,'XColor',[1 1 1]*0.7);
        subplot(2,2,2)
        imagesc(det); title('detected');
        set(gca,'Parent',hFig,'YColor',[1 1 1]*0.7,'XColor',[1 1 1]*0.7);
        subplot(2,2,3)
        imagesc(img_overlay); title('overlay');
        copyfile(fn_spotWin.name,['..\detectedSpots.mat' ]);
        set(gca,'Parent',hFig,'YColor',[1 1 1]*0.7,'XColor',[1 1 1]*0.7);
        cd('..\')
    end

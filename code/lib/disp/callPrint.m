for i = 1:4 % display and print all four cases 
    set(hPopupData,'Value',i);
    isTrace = i;
    updScat;

    if genOverlayImg
        figure(figImg)
        labelTrLen=[];
        if isTrace>=2
            labelTrLen = sprintf('_minTrLen%.02fms',minTrLen*acqTime*1000);
        end
        labelTime = sprintf('%.02fms_%.02fsec_1stFrm%i',acqTime*1000,acqTime*numFrames*binFrame,fr1);

        fOut = sprintf('recruitment%s%s_%s.tif',traceORspot{isTrace},labelTrLen,labelTime);
        imgOut = getframe(gcf);
        imwrite(imgOut.cdata,fOut,'Compression', 'none')
    end

    cfg.overlay.fr1 = fr1;
    cfg.overlay.fr2 = fr2;
    cfg.overlay.frames = frames;
    fOut_{1} = sprintf('recruitment%s_%s.tif',traceORspot{1},labelTime);
    fOut_{2} = sprintf('recruitment%s%s_%s.tif',traceORspot{2},labelTrLen,labelTime);
    fOut_{3} = sprintf('recruitment%s%s_%s.tif',traceORspot{3},labelTrLen,labelTime);
    fOut_{4} = sprintf('recruitment%s%s_%s.tif',traceORspot{4},labelTrLen,labelTime);
    [CB,~,~,hFig] = plotColorBar([1:numFrames]*acqTime,6);
    set(gcf,'PaperPositionMode','auto')
    colorbarPrint = sprintf('%s_colorbar.tiff',fOut);
    print(colorbarPrint,'-dtiff','-r80'); 
end
close(figImg)
close(hFig)
    

%stack name
currDir = cd;
posParanth= strfind(currDir,'\');
txtDate = currDir(posParanth(end-3)+1:posParanth(end-2)-1)
txtCell = currDir(posParanth(end-2)+1:posParanth(end-1)-1)
txtAcqNo = currDir(posParanth(end-1)+1:posParanth(end)-1)
txtCoeff = currDir(posParanth(end)+6:end)
foutStack = sprintf('overlayStack_%s_%s_%s_%s_%s_%s.tif',txtDate,txtCell,txtAcqNo,txtCoeff,labelTrLen,labelTime)

%% write to stack
overlayImg_1 = uint8(imread(fOut_{1}));
imwrite(overlayImg_1,foutStack,'Compression','none')
overlayImg_2 = uint8(imread(fOut_{2}));
imwrite(overlayImg_2,foutStack,'Compression','none','WriteMode','append')
overlayImg_3 = uint8(imread(fOut_{3}));
imwrite(overlayImg_3,foutStack,'Compression','none','WriteMode','append')
overlayImg_4 = uint8(imread(fOut_{4}));
imwrite(overlayImg_4,foutStack,'Compression','none','WriteMode','append')
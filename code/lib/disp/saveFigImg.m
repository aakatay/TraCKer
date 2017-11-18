%saves the current figure image to saveFigImgFN
imgFig = getframe(gcf); 
imgOut = imgFig.cdata;
imwrite(imgOut,saveFigImgFN,'WriteMode','append','Compression', 'none') 
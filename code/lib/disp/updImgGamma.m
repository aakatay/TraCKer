function updImgGamma(hObject,event,nTileXY, hTextImgGamma, hImg,img,imgMax_,hImgMin, hImgMax, hLineGamma)
    ydata=get(hLineGamma,'Ydata'); mxHist = ydata(end);
    hImgGamma = hObject;

    scaleImg;
    set(hImg,'Cdata',img);
    
    %gamma line
    nP = (imgMax-imgMin)/imgMax*1000;
    gammaLine = 1:nP;
    gamma = gammaLine.^imgGamma;
    gamma = gamma/max(gamma)*mxHist;
    set(hLineGamma,'XData',imgMin+gammaLine/1000*imgMax,'YData',gamma);
    set(hTextImgGamma,'String',sprintf('gamma: %.02f',imgGamma));
    

    
end

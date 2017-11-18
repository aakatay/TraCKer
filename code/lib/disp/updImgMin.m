function updImgMin(hObject,event,nTileXY,hTextImgMin, hImg,img,imgMax_,hImgMax, hImgGamma,hLineMin,hLineGamma)
    hImgMin = hObject;
    scaleImg;
    set(hLineMin,'XData',[imgMin imgMin]);
    set(hImg,'Cdata',img);
    set(hTextImgMin,'String',sprintf('min: %.02f',imgMin));
end

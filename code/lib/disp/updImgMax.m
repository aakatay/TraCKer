function updImgMax(hObject,event,nTileXY,hTextImgMax, hImg,img,imgMax_,hImgMin, hImgGamma,hLineMax,hLineGamma)
    hImgMax = hObject;
    scaleImg;
    set(hLineMax,'XData',[imgMax imgMax]);
    set(hImg,'Cdata',img);
    set(hTextImgMax,'String',sprintf('max: %.02f',imgMax));
end
function updImg(hObject,event,hTextImgMin, hTextImgMax, hTextImgGamma, hImg,img,imgMin,imgMax,imgGamma)
    slideVal = get(hObject,'Value');
    
    findPeaks; % find [x y] after coeff update
    set(hImg,'Xdata',x,'Ydata',y)
    set(hTextImgMin,'String',sprintf('min: %i',imgMin))
    set(hTextImgMax,'String',sprintf('max: %i',imgMax))
    set(hTextImgGamma,'String',sprintf('gamma: %i',imgGamma))
    
end

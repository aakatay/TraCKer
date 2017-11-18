function updSpotNum(hObject,event,hTextSpotNum,numSpots)

    slideVal = get(hObject,'Value');
    spotNum = round((numSpots-1)*slideVal)+1;
    set(hTextSpotNum,'String',sprintf('spot#: %i',spotNum));

end
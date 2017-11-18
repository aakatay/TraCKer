function updScrollState(hButton,hTextSpotNumState)

    textBtn = get(hButton,'String');
    if strcmp(textBtn,'start scroll')
        set(hTextSpotNumState,'String','scroll'); % textState 
        set(hButton,'String','stop scroll'); % textBtn 
    elseif strcmp(textBtn,'stop scroll')
        set(hTextSpotNumState,'String','select spot'); % textState 
        set(hButton,'String','start scroll'); % textBtn 
    elseif strcmp(textBtn,'set sections')
        set(hTextSpotNumState,'String','set cross section'); %  
        set(hButton,'String','select spot'); %  
    elseif strcmp(textBtn,'select spot')
        set(hTextSpotNumState,'String','select spot'); %  
        set(hButton,'String','set sections'); %  
    end
    
end
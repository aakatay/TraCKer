function updFrames(hObject, event,hScat,hPopupData,hTextFrames,hFramesSliderText1,hFramesSliderText2,hNumFramesSliderText,hFramesSlider1,hFramesSlider2,hNumFrames,hPlay,acqTime,frames)
%function updFrames(hObject, event,hScat,hPopupData,hTextFrames,hFramesSliderText1,hFramesSliderText2,hNumFramesSliderText,hFramesSlider1,hFramesSlider2,hNumFrames,hPlay,acqTime,X,Y,TraceX,TraceY,trInf2,frmNoSpot,frmNoTrace,CM,q,frames)
    % read sliders
    fr1 = round(get(hFramesSlider1,'Value')); % current Frame
    fr2 = round(get(hFramesSlider2,'Value')); % last Frame
    fr1_ = round(get(hNumFrames,'Value')); % previous current Frame
    numFrames = fr2 - fr1+1;
    if fr1_+numFrames-1 > frames, fr1_ = frames-numFrames+1;  end
    if hObject == hFramesSlider1
        set(hNumFrames,'Value',fr1);
    elseif hObject == hNumFrames
        fr1 = fr1_;
        fr2 = fr1 + numFrames-1;
        set(hFramesSlider1,'Value',fr1);
        set(hFramesSlider2,'Value',fr2);
    end   
    
    % time
    t=numFrames*acqTime;
    t1=fr1*acqTime;
    t2=t1+t;
    set(hTextFrames,'String',sprintf('time:%0.1f(%0.1f-%0.1f)',t,t1,t2));
    
    % update sliders
    set(hFramesSliderText1,'String',sprintf('fr1:%i',fr1));
    set(hFramesSliderText2,'String',sprintf('fr2:%i',fr2));
    set(hNumFramesSliderText,'String',sprintf('#fr:%i',numFrames));
    
    if hObject == hFramesSlider2 || hObject == hNumFrames
        set(hPlay,'String','play');
    else
        set(hPlay,'String','play window');
    end
    
    % frame process
    isTrace = get(hPopupData,'Value');    
    updScat;
    
    
end

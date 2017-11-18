function updBckCoeff(hObject,event,hImg,A,maxA,hTextBckCoeff)
    slideVal = get(hObject,'Value');
    Coeff = (slideVal*10)+1;
    
    Adisp = uint16(A-maxA/Coeff);
    set(hImg,'CData',Adisp)
    set(hTextBckCoeff,'String',sprintf('%.02f',Coeff));
    
    
end

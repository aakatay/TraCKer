function updXsection(hObject,event,hCrossSectX,hCrossSectY,hSect,IdatSHOW,Ifit2SHOW_2,BigWindowSize,intPeak)

    hSectDataX = hSect(1);
    hSectDataY = hSect(2);
    hSectFitX = hSect(3);
    hSectFitY = hSect(4);
    
    slideValX = get(hCrossSectX,'Value');
    sectPosY = slideValX*(BigWindowSize-1)+1;
    slideValY = get(hCrossSectY,'Value');
    sectPosX = slideValY*(BigWindowSize-1)+1;
    posSect = [sectPosX sectPosY sectPosX sectPosY]/2;
    if hCrossSectX == hObject % x section
        [sectDataX,sectDataY,sectFitX,sectFitY] = getCrossSect(IdatSHOW,Ifit2SHOW_2,posSect,intPeak)
        set(hSectDataX,'YData',sectDataX);
        set(hSectFitX,'YData',sectFitX);
    else % y section
        [sectDataX,sectDataY,sectFitX,sectFitY] = getCrossSect(IdatSHOW,Ifit2SHOW_2,posSect,intPeak)
        set(hSectDataY,'YData',sectDataY);
        set(hSectFitY,'YData',sectFitY);        
    end

end
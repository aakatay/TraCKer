function [mag, pos, m, n] = calcMaxMag(img,mag)
% calculates the maximum possible magnification of an image for display

    sizeNotOK = 1;    
    while sizeNotOK    
        v=size(img).*mag;
        m=v(2);
        n=v(1);
        pos_=get(0,'ScreenSize');
        pos = pos_(3:4); 
        check=pos - [m n+40];
        if isempty(find(check<0))
            sizeNotOK = 0;
        else
            mag= mag-1;
        end
       %mag= mag-1;
    end    
    pos = pos - [m n-40];
end
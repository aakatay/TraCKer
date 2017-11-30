function [mag, pos, m, n] = calcMaxMag(varargin)
% calculates the maximum possible magnification of an image for display
    if nargin == 2
        img = varargin{1};
        mag = varargin{2};
    elseif nargin == 3
        img = varargin{1};
        mag = varargin{2};
        scrnSizeIn = varargin{3};
    end

    sizeNotOK = 1;    
    while sizeNotOK    
        v=size(img).*mag;
        m=v(2);
        n=v(1);
        if exist('scrnSizeIn')
            pos=scrnSizeIn;
        else
            pos_=get(0,'ScreenSize');
            pos = pos_(3:4); 
        end
        check=pos - [m n+40];
        if isempty(find(check<0))
            sizeNotOK = 0;
        else
            mag= mag-1;
        end
       %mag= mag-1;
    end    
    pos_=get(0,'ScreenSize');
    pos = pos_(3:4); 
    pos = pos - [m+10 n+20];
end
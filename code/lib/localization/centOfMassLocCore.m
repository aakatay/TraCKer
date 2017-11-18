function [Xc,Yc] = centOfMassLocCore(Window,WindowSize)
%FIND Intensity Center
    TopX=0;
    TopY=0;
    TopColum=0;
    TopRow=0;
    WSumX=0;
    WSumY=0;

    for j=1:WindowSize
       TopX(j)=sum(Window(:,j));
    end
    TopX=TopX-min(TopX);
    TopRow=sum(TopX);

    for j=1:WindowSize
        WSumX=WSumX+j*TopX(j);
    end

    for ii=1:WindowSize
       TopY(ii)=sum(Window(ii,:));
    end
    TopY=TopY-min(TopY);
    TopColum=sum(TopY);

    for ii=1:WindowSize
        WSumY=WSumY+ii*TopY(ii);
    end

    Xc=WSumX/TopRow;
    Yc=WSumY/TopColum;
end  
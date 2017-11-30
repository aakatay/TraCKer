% center of mass localization
% B --> (X_, Y_, INT_)
% B: boundaries of the BW peaks
% (X_, Y_): coordinates of the localized peak
outsideCall=0;
if exist('B')
    c=cell2mat(B(m));
    Py=uint16(mean(c(:,1)));
    Px=uint16(mean(c(:,2)));
else %  Px Py defined outside
    outsideCall=1; 
    k=1;
end
Px0 = double(Px);
Py0 = double(Py);

% adjust the big window center position for the spots at
% the edges
PxBW = Px;
PyBW = Py;

if (Px-(BigWindowSize+1)/2)<1
    PxBW=(BigWindowSize+1)/2;
end
if (Py-(BigWindowSize+1)/2)<1
    PyBW=(BigWindowSize+1)/2;
end
if (Px+(BigWindowSize+1)/2)>En1*div
    PxBW=En1*div-(BigWindowSize+1)/2;
end
if (Py+(BigWindowSize+1)/2)>Boy1*div
    PyBW=Boy1*div-(BigWindowSize+1)/2;
end
dPy = double(PyBW)-double(Py);
dPx = double(PxBW)-double(Px);
Px = PxBW; Py = PyBW;
%DEFINE Window
yy = Py-(WindowSize+1)/2;
y1 = yy + 1; y2 = yy + WindowSize;
xx = Px-(WindowSize+1)/2;
x1 = xx + 1; x2 = xx + WindowSize;
Window=IMG(y1:y2,x1:x2);
%if isdbg,figure(890); imagesc(Window);axis image; end

%DEFINE Big Window
yy = PyBW-(BigWindowSize+1)/2;
y1 = yy + 1; y2 = yy + BigWindowSize;
xx = PxBW-(BigWindowSize+1)/2;
x1 = xx + 1; x2 = xx + BigWindowSize;
BigWindow=IMG(y1:y2,x1:x2);

%background
BACKmean=[min(mean(BigWindow,1)),min(mean(BigWindow,2))];
BACK =min(BACKmean);

%FIND Total Intensity
INT_=sum(sum(Window))-BACK * (WindowSize)^2;


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

Xc(k)=WSumX/TopRow;
Yc(k)=WSumY/TopColum;
if isnan(Xc(k)), Xc(k)=(WindowSize+1)/2; end;
if isnan(Yc(k)), Yc(k)=(WindowSize+1)/2; end;

PXc=uint8(Xc(k));
PYc=uint8(Yc(k));

%center of intensity
X_=double(Px)+Xc(k)-double((WindowSize+1)/2);
Y_=double(Py)+Yc(k)-double((WindowSize+1)/2);

X_ = X_-dPx;
Y_ = Y_-dPy;


%ddxx = [ddxx double(Px0) - X_];
%ddyy = [ddyy double(Py0) - Y_];
ccc=4;
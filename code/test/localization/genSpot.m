clear
fname = 'a.tif';
save('fname','fname')

noiseRat = 0.5;
SNR = 1/noiseRat;
intMax = 5000;

dx = 1; dy = -1;
Y = [2  20  10+dy  18  3];
X = [10  5  10+dx  14  14];
DY = rand(size(X))-0.5;
DX = rand(size(X))-0.5;

Y2 = Y + DY;
X2 = X + DX;

%X = X-4; Y = Y-4;
P = [1 1 1 1 1]; % intensity
N = 32;
a = zeros(N,N);
gaus2 = fspecial('gaus',9,0.7);
gaus2 = gaus2/max(gaus2(:));
for i = 1:numel(P)
    y=Y(i);x=X(i);p=P(i);
    dy=DY(i);dx=DX(i);p=P(i);
    %gausSpot = p*exp(-((dx).^2 ...
    %                  +(dy).^2)/s1^2/2);
    a(y:y+8,x:x+8) = gaus2*p+ a(y:y+8,x:x+8);
end

%save('XY','Xgen','Ygen');


mx = max(a(:));
mnP = min(P);

gausBlur = fspecial('gaus',3,2);
noise = rand(size(a));
noise = imfilter(noise,gausBlur,'symmetric');
mxNoise = max(noise(:));
noise = noise*mnP/SNR/mxNoise;
a = noise + a ;
imagesc(a)
axis equal, axis tight
imwrite(uint16(a*intMax),'a.tif');


fun = @(c,x) c(1)+c(2)*exp(-((x(:,1)-c(3))/c(4)/sqrt(2)).^2-((x(:,2)-c(5))/c(6)/sqrt(2)).^2);
c0 = [0 1 5 2 5 2];
sp = a(6:14,6:14);
[cc1,err1G5by5,err1res_,EXITFLAG]=gausFit(sp,fun,c0);
cc1
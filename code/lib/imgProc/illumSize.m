fn='AVG_illum_cell12_frm1-35.tif';
fn2 = 'illumBinary.tif';
fn3 = 'illumMass.tif';
th1 = 185;
th2 = 450;
illumPower = 1; % [mw]
P_D63=0.632;
radTh = 1-1/exp(1); % D63 diameter
illumPower = illumPower*P_D63;
%radTh = 0.5; % FWHM

TH = (th1+(th2-th1)*radTh);




A=double(imread(fn));
Amx = max(A(:));

Ab = im2bw(A/Amx,TH/Amx);

[B,L] = bwboundaries(Ab,'noholes');

for i = 1: numel(B)
    sz(i)=numel([B{i}]);
end
[v,ix]=max(sz);

AL = L;
AL(AL~=ix)=0;
imagesc(AL)
L = watershed(AL);
L=im2bw(L-1,0);
%L(L > 0)=1;


figure(2);
imagesc(L)
imwrite(Ab,fn2,'Compression', 'none'); 
%%
figure(3);
L2 = conv2(double(L),ones(3),'same');
L2 = im2bw(L2,0);
imagesc(L2)

imwrite(L2,fn3,'Compression', 'none'); 
IllumSz = sum(L2(:));
%IllumSz = 512^2;
%%
pxSz = 160/1.5*1e-3; % um
pxArea = pxSz^2;
illumMicroMeterSqr = IllumSz*pxArea;
illumDensity = illumPower/illumMicroMeterSqr*1000; % [uw/um^2]

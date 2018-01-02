% called by callFluorescenceSpectra.m
% displays intensity bleed throughs
close all;
load('channelIntensity'); % 'chls','chfl','flx'

% output
imgFout = 'SNRchannelLaserFluo.tif';
MAToutFN  = 'SNRchannelLaserFluo.mat';
if isBRT
    imgFout = ['BRT-' imgFout];
    MAToutFN = ['BRT-' MAToutFN];
else
    imgFout = ['PRB-' imgFout];
    MAToutFN = ['PRB-' MAToutFN];
end
   
if isCy3
    imgFout = ['CY3-' imgFout];
    MAToutFN = ['CY3-' MAToutFN];
else
    imgFout = ['GFP-' imgFout];
    MAToutFN = ['GFP-' MAToutFN];
end 

% no 4th laser and 4th fluo
chls4by4 = nan(4);
chfl4by4 = nan(4);
flx4by4  = nan(4);
if isBRT
    nCH = 4;
    chls4by4(1:3,1:3) = chls;
    chfl4by4(1:3,1:3) = chfl;
else
    nCH = 1;
    chls4by4(1,1:3) = chls;
    chfl4by4(1,1:3) = chfl;
end
flx4by4(1:3,1:3)  = flx;
flB = fl.B; % brightness
lsP = ls.P; % laser power
flB(4) = 0;
lsP(4) = 0;

% 64 square grid images
w = 32; % width
T0 = repmat([ones(w,1) zeros(w,1)],1,w/2);
T0(:,w/2+1)=0;
T0(1,:)=-1;
T0(:,end)=-1;
%figure(61);imagesc(T0)

a=[4 3; 1 2];
a=[a+12 a+8;a a+4];
a=[a+48 a+32;a a+16];
A=repelem(a,w,w);
A = A*1i;
am = a*1i;
ax = a*1i;

% 
G=imread('spectrumData\channelLaserFluoGRID.tif');
figure(62);imagesc(G); axis image;


for i = 1:nCH % channels
    fluoEff = sys.fluoCollect; % fluo imaging efficiency
    ls = chls4by4(i,:).*lsP; % laser source leak
    fl = chfl4by4(i,:).*flB*fluoEff; % fluo
    FL = repmat(fl,4,1);
    t = flx4by4.*FL;
    for j = 1:4 % lasers
        for k = 1:4 % fluos
            T = T0;
            p = k+(j-1)*4+(i-1)*16;
            T(T0==0)=ls(j);T(T0==1)=t(j,k);
            %T(T0==0)=p;T(T0==1)=p;
            A(A==p*1i)=T;
            am(am==p*1i) = t(j,k);
            ax(ax==p*1i) = ls(j);
        end
    end
    dbg=0;
    if dbg
        Adbg = A;
        Adbg(isnan(A)) = 0;
        Adbg(A==-1) = 0;
        figure(63);imagesc(real(Adbg));
        colormap('gray')
        cc=1;
    end
end
%% display
% each combination (leaks)
A = real(A);
ax = real(ax);
am = real(am);
figure(63)
mxA = max(A(:));
A(A==-1)=mxA;
A(128:130,:)=mxA/1.5;
A(:,127:129)=mxA/1.5;
imagesc(A); axis image
colormap('gray')

%% save 
        
if exist(imgFout),delete(imgFout); end
figure(63)
imgFig = getframe(gcf);
imgOut = imgFig.cdata;
imwrite(imgOut,imgFout,'WriteMode','append','Compression', 'none') 
figure(62)
imgFig = getframe(gcf);
imgOut = imgFig.cdata;
imwrite(imgOut,imgFout,'WriteMode','append','Compression', 'none') 

save(MAToutFN,'am','ax','A');
if 0
    % sort flx
    a=[4 3; 1 2];
    a=[a+12 a+8;a a+4];
    [~,ix]=sort(a(:));
    flxS(ix)=flx(:);
end
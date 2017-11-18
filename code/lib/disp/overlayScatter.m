close all;
clear

WindowSize = 5;
nWalkAv = 3;
pxSz = 100; % [nm] pixel size
recRes = 25; % [nm] recruitment spacing (forced resolution)


% 70*40*180 *4*4 (uint8)=>~80MB (uint16)=>~200MB (double)=>~720MB 
%% load files
imgLongExposure=dir('..\imgs\img*.tif');
fn = imgLongExposure.name;
a=imread(sprintf('..\\imgs\\%s',fn));
cropSz = (WindowSize+1)/2;
a = a(cropSz+1 : end-cropSz,cropSz+1 : end-cropSz  );
data = dir('xyzDataGausFilt*.mat');
load(data.name);
frames = size(X,2);
% load img file
if exist('fname.mat')
    load fname
    if ~exist(fname)
        if exist(sprintf('..\\%s',fname))
            fname = sprintf('..\\%s',fname);
        end
    end
end
imgFin = fname;
if ~exist(imgFin), imgFin2 = rdir(sprintf('..\\%s',imgFin)); imgFin=imgFin2.name; end;
I = imread(imgFin);
[yNpx,xNpx]=size(I);
imgMax= max(I(:));

%% calculate maximum magnification
mag = pxSz/recRes;
mag = 5;
gSz = 5*mag;
gSig = mag;
gaus=fspecial('gaussian', gSz, gSig); % convolution function
mxGaus = max(gaus(:));
[mag, pos, m, n ] = calcMaxMag(a,mag);

%% recruitment detection
isVolDetect  = 1;
hWB = waitbar(0,'detecting spots');
if isVolDetect
    Awalk = uint8(zeros(mag*yNpx,mag*xNpx,nWalkAv));
    w = 1;
    %nf = 20; 
    nf= 1;
    for fr = 1:frames/nf % # of frames
        A = uint8(zeros(mag*yNpx,mag*xNpx));
        
        for tr = 1:sum(X(:,fr)~=0)
            x = X(tr,fr);
            y = Y(tr,fr);
            x = ceil(mag*x);
            y = ceil(mag*y);
            A(y,x) = 1; % pits
            %imagesc(A)
            %pause(0.2)
        end
        Awalk(:,:,1:nWalkAv-1) = Awalk(:,:,2:nWalkAv);
        Awalk(:,:,nWalkAv) = A;
        if fr>=nWalkAv
            AwalkAv(:,:,w) = mean(Awalk,3); 
            convVol(:,:,w) = imfilter(AwalkAv(:,:,w),gaus,'symmetric');
            w = w + 1; 
        end
        waitbar(fr/frames*nf);
    end
    close(hWB);
    %  normalize data
    x1 = 1; y1 = 1;
    convVol = convVol(x1:end,y1:end,:);

    convVol = convVol-min(convVol(:));
    maxData = max(convVol(:));
    imgMul = (2^16-1)/maxData;
    mxGausIm = mxGaus*imgMul;
    convVol = convVol*imgMul;
    %stackWrite(convVol,'convVol.tif');
    %stackWrite(AwalkAv,'AwalkAv.tif');
    colormap gray;    imagesc(sum(convVol,3)); axis image
end
    %save('3Ddata','convVol','AwalkAv');

%% scatter overlay
tit = 'reconstructed image';
CM = genColorMap('jet',frames);
fig=figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos/2 m n]);
%axis image
axe=axes('Parent',fig,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off');
imagesc(a);
hold on;
X = X - cropSz;
Y = Y - cropSz;

isColorTime = 1;
if isColorTime
    for i = 1: frames
        CL = repmat(CM(i,:),size(X,1),1);
        scatter(X(:,i),Y(:,i),1,CL,'.');
        %pause
    end
else
    scatter(X(:),Y(:),1,'.');
end
hold off;
set(gcf,'PaperPositionMode','auto')

resLow = round(1000/mag/20)*20;
resHigh = round(1000/mag/20)*20*2;
fnLow = sprintf('recruitment_mag%i_R%i.tiff',mag,resLow);
fnHigh = sprintf('recruitment_mag%i_R%i.tiff',mag,resHigh);
rLow = sprintf('-r%i',resLow);
rHigh = sprintf('-r%i',resHigh);
%print('recruitment.bmp','-dbmp16m'); 
%print('recruitment.bmp','-dbmp16m'); 
%print('recruitment.eps','-deps'); 
%print('recruitment.emf','-dmeta'); 
print(fnLow,'-dtiff',rLow); 
print(fnHigh,'-dtiff',rHigh); 
	
%print('picture1','-dtiff','-r300')
%print('-dpsc2','-zbuffer','-r200')


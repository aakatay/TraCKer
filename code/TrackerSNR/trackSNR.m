% reads the 'filename_001.tif' and data from tracking of  'filename_002.tif' 
% processes 'filename_001.tif'

clear all;
close all;

%% input files
load fname;
fname0 = [fname(1:end-5) '1.tif'];
tracesFN = rdir('*\traceData0-coeff*');
traces2FN = rdir('*\traceJmplessData-coeff*');

%% output files
SNRmovieFN = 'SNRmovie.tif';
SNRmovieConvFN = 'SNRmovieConv.tif';
SNRplotFN = 'SNRplot.tif';
numSMplotFN = 'SNRnumSMplot.tif';
snrSTDplotFN = 'SNRstdplot.tif';
SMdataFN = 'SMdata.mat';

%%
isRTdisp = 0;
stdWin = 10; % number of frames to calc. std
SNRmul = 1000; % for 'SNRimg.tif'
acqTime = 0.2;
fcall = 'trackSNR';

%% read acquisition files
imgInf = imfinfo(fname);
nfr = numel(imgInf);
imgInf = imfinfo(fname0);
nfr0 = numel(imgInf);
waWinSz = nfr0-nfr+1; % walking average window size
szXY = [imgInf(1).Height imgInf(1).Width];
xc1 = 160;
yc1 = szXY(2)-110;
xc2 = xc1+20;
yc2 = yc1+25;

xc1 = 1;
yc1 = szXY(2)-1;
xc2 = xc1+100;
yc2 = yc1+100;

xc1 = 1;
yc1 = 1;
xc2 = xc1+50;
yc2 = yc1+50;


%% read trace data
load(tracesFN.name); % trInf 
load(traces2FN.name); % TraceX2 TraceY2
X = TraceX2;
Y = TraceY2;
%trInf = trInf(515:518,:);

frm1 = trInf(:,1);
%trInf(frm1~=1)=[]; % remove new recruitments

n = size(trInf,1); % # of single molecule
frm1 = trInf(:,1);
frm2 = trInf(:,2)+trInf(:,1)-1;
xyIx = trInf(:,3); % index to X Y arrays
minNumFrm = nfr; % single molecule selection
wsz = 6; % window size
s = wsz/2; 

tit = 'SNR image';
mag = 1; % display size
[mag, pos, szx, szy ] = calcMaxMag(zeros(szXY),mag);
szXYmag = [szx, szy];
colormap('gray');
pos(1) = pos(1) - 800;
pos(1) = pos(1)- 700;
if isRTdisp
    figSNR = figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos/2 szXYmag(1) szXYmag(2)]);
    axeSNR = axes('Parent',figSNR,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren','XLim',0.5+[0 szXYmag(1)],'YLim',0.5+[0 szXYmag(2)]);
end

dbg = 0;
if dbg    
    dbgImg = zeros(szXYmag);
    xdbg = round(trInf(:,4));
    ydbg = round(trInf(:,5));
    dbgImg(sub2ind(size(dbgImg),ydbg,xdbg)) = 1;
    figure(209); imagesc(dbgImg); axis image
    imwrite(dbgImg,'dbgImg.tif')
end

%% crop single molecule windows
f = 1; % a index
aIX = cell(n,1);
fIX = cell(n,1);
snrImg = zeros(szXYmag);
m = 1; % SM index for SNR
tt =0;
for i = 1:nfr % each frame
    snrIMG = zeros(szXYmag);
    A = imread(fname0,i);
    A = padarray(A,[s s]);
    IX = find((frm1<=i) .* (i<=frm2));
    fTrace = i - frm1(IX) + 1;
    
    x = ceil(X(xyIx(IX)+fTrace-1));
    y = ceil(Y(xyIx(IX)+fTrace-1));    
    
    
    % find non-overlapping SM
    smMap = zeros(szXY);
    smMap(sub2ind(szXY,y,x)) = 1;
    smMapCv = conv2(smMap,ones(wsz),'same');
    smMapCv(smMapCv~=1) = -10000;
    smMapCv2 = conv2(smMapCv,ones(wsz),'same');
    smMapCv3 = zeros(szXY);
    smMapCv3(smMapCv2 == wsz^2) = 1;
    smMapCv4 = circshift(smMapCv3,[1 1]);
    smMapCv5 = smMapCv4+smMap;
    ixImg = find(smMapCv5==2);
    ixIMG = sub2ind(szXY,y,x)';
    
    dbg = 0;
    if dbg 
        figure(344);imagesc(smMap)
        figure(345);imagesc(smMapCv)
        figure(346);imagesc(smMapCv2)
        figure(347);imagesc(smMapCv3)
        figure(348);imagesc(smMapCv4)
        figure(349);imagesc(smMapCv5)
    end
    num_ixImg(i) = numel(ixImg);
    if isempty(ixImg), continue; end;
    
    ixSM_ = find(ismember(ixIMG,ixImg));
    ixSM = IX(ixSM_);
    
    for j = 1:numel(ixSM) % each single molecule in that frame
        
        ixsm = ixSM(j); % single molecule index
        y_ = y(j);
        x_ = x(j);
        if ((x_<xc1) || (y_<yc1) || (x_>xc2) || (y_>yc2)),
            tt = tt + 1;
            %continue; 
        end;
        if ((x_<=s) || (y_<=s) || (x_+s-1>szXY(1)) || (y_+s-1>szXY(2))), continue; end;
        a(:,:,f) = A(y_-s:y_+s-1,x_-s:x_+s-1); % crop
        aIX{ixsm} = [aIX{ixsm} f]; % index for the frame of the single molecule
        fIX{ixsm} = [fIX{ixsm} i]; % frames
        f=f+1;
    end
    
    if i >= stdWin % calculate SNR
        aIXix_ = find(~cellfun(@isempty,aIX)); % all SM found so far
        aIXix = aIXix_;
        aIXix( (   (i-frm1(aIXix)+1)<stdWin   )) = []; % remove SM shorter than stdWin
        aIXix( (   frm2(aIXix) < i  )) = []; % remove SM bleach before the current frame (i)
        mlast = m;
        snrImg = zeros(szXYmag);
        num_ixImg2(i) = 0;
        for j = 1:numel(aIXix) % for each SM 
            ix = aIXix(j);
            fr = fIX{ix}; % frames
            if fr(end) ~= i, continue; end; % no data in this frame
            if numel(fr)-stdWin+1 < 1, continue; end; % not enough data
            if fr(end)-stdWin+1 ~= fr(end-stdWin+1) % missing data in the last stdWin frames
                continue;
            end
            num_ixImg2(i) = num_ixImg2(i) + 1;
            aix = aIX{ix};
            asm = a(:,:,aix); % SM crop stack
            smLast = asm(:,:,end); % last image
            
            intWin = double(smLast(s-1:s+2,s-1:s+2));
            int = sum(intWin(:));
            peak = max(intWin(:));
            asm(s-1:s+2,s-1:s+2,:) = nan;
            a2D = double(reshape(asm,wsz^2,size(asm,3)));
             
            astd = std(a2D,0,1);
            astd(isnan(astd)) = [];
            b = mean(astd); % background std
            snr = peak/sqrt(peak+b^2);
            
            fTrace = i - frm1(ix) + 1;
            xyix = trInf(ix,3)+fTrace; % index to X Y arrays 
            Xs(m) = X(xyix);
            Ys(m) = Y(xyix);
            Frm(m) = i;
            TRinf(m,:) = trInf(ix,:);     
            B(m) = b;
            INT(m) = int;
            SNR(m) = snr;                   
            MIX(m) = ix; % SM index
            
            Xmag = round(Xs(m)*mag);
            Ymag = round(Ys(m)*mag);
            
            snrImg(Ymag,Xmag,j) = SNR(m);
            
            m = m + 1;
            %BCK(i) = mean(bckgrnd(:));
            %INT(i) = max(peak(:));
            %RECT(i,:) = rect;
        end
        snrIMG = sum(snrImg,3);
        sn{i} = find(snrImg~=0);
    end
    if isRTdisp
        figure(figSNR);
        set(axeSNR,'Parent',figSNR,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','XLim',1+[0 szXYmag(1)],'YLim',1+[0 szXYmag(2)]);
        imagesc(flipud(snrIMG))
        pause(acqTime)
    end
    SNRimg(:,:,i) = snrIMG;
    
end
%% save
save(SMdataFN,'Xs','Ys','Frm','TRinf','B','INT','SNR','MIX','aIX','fIX');
%sum(reshape(SNRimg,260^2,226),1)
if ~isempty(find(SNRimg*SNRmul>=2^16)), warning('SNR image saturated use lower a SNRmul');end

%% SNR movie
SNRmov = SNRimg*SNRmul;
stackWrite(SNRmov,SNRmovieFN);
if 0 
    cvWin = ones(5);
    SNRmovCv = convn(SNRmov,cvWin,'same');
    stackWrite(SNRmovCv,SNRmovieConvFN);
end

%% SNR plot
figure(201);
SNR2d = reshape(SNRimg,[szXYmag(2)*szXYmag(1) nfr]);
SNR2dSUM = sum(SNR2d,1);
SNR2dNUM = sum(SNR2d>0,1);
snrMean = SNR2dSUM./SNR2dNUM;

plot(snrMean); 
title('SNR mean')
xlabel('frames')
ylabel('SNR')   
imgFig = getframe(gcf);
imgOut = imgFig.cdata;
imwrite(imgOut,SNRplotFN);

%% SNR std plot
figure(202);
SNR2d = reshape(SNRimg,[szXYmag(2)*szXYmag(1) nfr]);
SNR2dSTD = std(SNR2d,1);

plot(SNR2dSTD); 
title('SNR STD')
xlabel('frames')
ylabel('SNR STD')   
imgFig = getframe(gcf);
imgOut = imgFig.cdata;
imwrite(imgOut,snrSTDplotFN);

%% number of SM plot
figure(204);
plot(num_ixImg2);
title('number of SM for SNR reading')
xlabel('frames')
ylabel('Number of Single Molecules')   
imgFig = getframe(gcf);
imgOut = imgFig.cdata;
imwrite(imgOut,numSMplotFN);


return;

%% intensity traces
figure(210); maximize;clf;
hold on;
%load(SMdataFN); % 'Xs','Ys','Frm','TRinf','B','INT','SNR','MIX');
smix = unique(MIX); % # selected SM
if smix(1)==0, smix(1)=[]; end; % remove zero

smixSel = [11 14 16 21];
ns = numel(smix);
ns2 = 0;
ns3 = 0;
FRM = [];
intShft = 10e4;
intShft = 0;
ixs = [];
for i = 1:ns
    ix = smix(i); % SM index
    mix = find(MIX==ix);
    frm = Frm(mix);
    FRM(i)=numel(frm);
    if numel(frm)<100,continue; end;
    if numel(frm)<217,continue; end;
    
    ns2 = ns2 + 1;
    %if ~ismember(ns2,smixSel),continue; end;
    ns3= ns3+1;
    ixs = [ixs ix];
    int = INT(mix);
    INTsave(:,ns3) = int;
    subplot(4,1,ns3)
    plot(frm,int+intShft*(ns3-1));
    ylim([0 max(int)*1.2])
    pause(0.5)
end
save('SNRintTraces','INTsave','frm','ixs')
savefig('SNRintTraces');
hold off;

if dbg
    ix = find(~cellfun(@isempty,aIX));
    ix(1:numel(MIX),2)=MIX';
    
end

if 0  % write SM stack
    %%
    ix = find(~cellfun(@isempty,aIX));
    ixs=ix(1);
    frm = aIX{ixs};
    a(:,:,frm);
    xs = round(trInf(ixs,4));
    ys = round(trInf(ixs,5));
    fname = sprintf('SM%03i_x%iy%i',ixs,xs,ys);
    stackWrite(a,fname)
    
end

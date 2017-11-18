% finds distance based colocalization of two localization sets
clear all;
xyzDIR  = rdir('*\**\xyzDataGaus-coeff0*.mat');
trInfDIR  = rdir('*\**\traceData0-coeff0*.mat');

isDbg =0;
isTrackData =1; % select data type
fDbg = [1:254];

mag = 3;
if isTrackData
    load(trInfDIR(1).name); % CLC
    tr1 = trInf;
    load(trInfDIR(2).name); % DYN
    tr2 = trInf;
    
    % CLC
    x1 = tr1(:,4);
    y1 = tr1(:,5);   
    X1 = floor(x1*mag+1); Y1 = floor(y1*mag+1); 
    INT1 = tr1(:,6);
    F1 = tr1(:,1) + tr1(:,2) - 1;
    
    % DYN
    x2 = tr2(:,4);
    y2 = tr2(:,5);
    X2 = floor(x2*mag+1); Y2 = floor(y2*mag+1); 
    INT2 = tr2(:,6);
    F2 = tr2(:,1) + tr2(:,2) - 1;
else % localization data
    load(xyzDIR(1).name); % CLC
    F1 = frmNoSpot;
    if isDbg, X = X(find(ismember(F1,fDbg))); Y = Y(find(ismember(F1,fDbg)));F1=F1(find(ismember(F1,fDbg)))-fDbg(1)+1; end;
    x1 = X; y1 = Y;
    INT1= INT; X1 = floor(X*mag+1); Y1 = floor(Y*mag+1); 
    load(xyzDIR(2).name); % DYN
    F2 = frmNoSpot;
    if isDbg, X = X(find(ismember(F2,fDbg))); Y = Y(find(ismember(F2,fDbg)));F2=F2(find(ismember(F2,fDbg)))-fDbg(1)+1; end
    x2 = X; y2 = Y;
    INT2= INT; X2 = floor(X*mag+1); Y2 = floor(Y*mag+1); 
end

szY = cfg.img.boy*mag;
szX = cfg.img.en*mag;
szT = cfg.img.lastFrm;
if isDbg, szT = numel(fDbg); end
sz = [szY szX szT];

ix1 = sub2ind(sz,Y1,X1,F1);
ix2 = sub2ind(sz,Y2,X2,F2);
% sort
[~,ix1s] = sort(ix1);
[~,ix2s] = sort(ix2);
X1=X1(ix1s);
Y1=Y1(ix1s);
F1=F1(ix1s);
INT1=INT1(ix1s);
x1=x1(ix1s);
y1=y1(ix1s);

X2=X2(ix2s);
Y2=Y2(ix2s);
F2=F2(ix2s);
INT2=INT2(ix2s);
x2=x2(ix2s);
y2=y2(ix2s);



%C(sub2ind(sz,Y1,X1,F1)) = INT1;
%D(sub2ind(sz,Y2,X2,F2)) = INT2;

fname1 = 'diffLocCLC.tif';
fname2 = 'diffLocDYN.tif';
fname12 = 'diffLocCLC-DYN.tif';
%stackWrite(C,fname1);
%stackWrite(D,fname2);

%delete(fname12);
C = zeros(szY,szX,1);
D = zeros(szY,szX,1);
Dprev = 0;
DixPrev=zeros(size(D,1),size(D,2));
ixs = 10000; % index shift (to disting. same channel overlaps)
DnzNlast=ixs;
CnzNlast=ixs;
for j = 1:szT
    [ixs1,~] = find(F1==(j-1));
    [ixs2,~] = find(F2==j);
    C = zeros(szY,szX,1);
    D = zeros(szY,szX,1);
    C(sub2ind(sz,Y1(ixs1),X1(ixs1),repmat(1,size(ixs1)))) = 1;
    D(sub2ind(sz,Y2(ixs2),X2(ixs2),repmat(1,size(ixs2)))) = 1;
    clear imgCol
    Dix = D;
    Dix0 = Dix;
    [~,Dnz] = find(Dix(:)'>0); % non zero indices
    DnzN = numel(Dnz);
    Dix = DixPrev;
    Dix(Dnz) = DnzNlast+[1:DnzN]; 
    Dix0(Dnz) = DnzNlast+[1:DnzN]; 
    DnzNlast = DnzNlast+DnzN;
    if j>1, Dprev = D0; DixPrev = Dix0; end;

    Cix = C;
    [~,Cnz] = find(Cix(:)'>0); % non zero indices
    CnzN = numel(Cnz);
    Cix(Cnz) = CnzNlast+[1:CnzN];
    CnzNlast = CnzNlast+CnzN;

    imgCol(:,:,1)=im2bw(D+Dprev); % red
    imgCol(:,:,2)=im2bw(C); % green
    imgCol(:,:,3)=0;
    %imwrite(uint16(imgCol),fname12,'WriteMode','append','Compression', 'none'); %
    imgCom(:,:,j) = Dix*1i + Cix;
    D0 = D;
    cccc=3;
end
cwszR = 1.5; % [px]
cwsz = ceil(cwszR*mag); % convolution window size
imgComConv = convn(imgCom,ones(cwsz),'same');
szimg = size(imgComConv);
cLix = find(conj(imgComConv)~=imgComConv&imag(imgComConv)*1i~=imgComConv); % close localizations
cL = imgComConv(cLix);
imgComConv2 = zeros(szimg);
imgComConv2(cLix) = 1;
%imgComConv = reshape(imgComConv,szimg);
[cL2,ixu,~] = unique(cL);
%[~,ixus] = sort(ixu); cL2 = cL2(ixus);

%cL2 = cL;
%%
cix = real(cL2)-ixs;
dix = imag(cL2)-ixs;
if 1
    ixdel = find(cix>numel(x1)|dix>numel(x2)|cix<1|dix<1);
    cix(ixdel)=[];
    dix(ixdel)=[];
end

xc = x1(cix); % clathrin
yc = y1(cix); 
xd = x2(dix); % dynamin
yd = y2(dix); 
d = sqrt((xc-xd).^2+(yc-yd).^2);


ncd=numel(cix);
figure(13); 
plotyy(1:ncd,cix',1:ncd,d)

hold on
plot(dix','r')
hold off;

% filter out large dist
dmax = 1.75;
ixIgnore = find(d>dmax);
nIgnore = numel(ixIgnore);
sprintf('%.1f%% is ignored',nIgnore/numel(d)*100)
d(ixIgnore) = [];
cix(ixIgnore) = [];
dix(ixIgnore) = [];


ncd=numel(cix);
figure(14); 
plotyy(1:ncd,cix',1:ncd,d)

hold on
plot(dix','r')
hold off;

xc = x1(cix); % clathrin
yc = y1(cix); 
xd = x2(dix); % dynamin
yd = y2(dix); 
fc = F1(cix);
fd = F2(dix);

figure(15); 
plotyy(1:ncd,fc,1:ncd,fd)

save('tandemTrack','xc','yc','xd','yd','fc','fd');



% reads the 'filename_001.tif' and 'filename_002.tif' 
% processes 'filename_001.tif'

clear all;
close all;
load fname;

%% read files
imgInf = imfinfo(fname);
nfr2 = numel(imgInf);
imgInf = imfinfo(fname);
nfr = numel(imgInf);
waWinSz = nfr-nfr2+1; % walking average window size
frames = waWinSz/2+1:n

fname = [fname(1:end-5) '1.tif'];
tracesFN = rdir('*\traceData0-coeff*');
traces2FN = rdir('*\traceJmplessData-coeff*');

%xyFN = rdir('*\xyzDataGaus-coeff*');
load(tracesFN.name); % trInf 
load(traces2FN.name); % TraceX2 TraceY2
%load(xyFN.name); % X Y INT
X = TraceX2;
Y = TraceY2;

n = size(trInf,1); % # of single molecule
frm1 = trInf(:,1);
frm2 = trInf(:,2)+trInf(:,1)-1;
xyIx = trInf(:,3); % index to X Y arrays
minNumFrm = nfr; % single molecule selection
wsz = 10; % window size
s = wsz/2; 

%% crop single molecule windows
f = 1; % a index
aIX = cell(n,1);
for i = 1:nfr % each frame
    A = imread(fname,i);
    A = padarray(A,[s s]);
    IX = find((frm1<=i) .* (i<=frm2));
    for j = 1:numel(IX) % each single molecule in that frame
        ix = IX(j); % single molecule index
        fTrace = i - frm1(ix) + 1;
        aIX{ix} = [aIX{ix} f]; % index for the frame of the single molecule
        x = ceil(X(xyIx(ix)+fTrace-1)+s);
        y = ceil(Y(xyIx(ix)+fTrace-1)+s);
        
        a(:,:,f) = A(y-s:y+s-1,x-s:x+s-1);
        f=f+1;
    end
end

%% averages of each single molecule
clear A;
FRM = [];
j=1;
ixSel =[]; % selected SM
for i = 1:n % each single molecule
    frm = aIX{i};
    if numel(frm)<minNumFrm, continue; end;
    FRM = [FRM frm 0];
    ixSel = [ixSel i];
    A(:,:,j) = mean(a(:,:,frm),3);
    j = j + 1;
end
stackWrite(A,'singleMoleculeMeanStack.tif');


%% write trace images
% add a black frame between single mol trace images
if 0
    a(:,:,end+1) = 0;
    FRM(FRM == 0)=size(a,3);

    asort = a(:,:,FRM);
    stackWrite(asort,'singleMoleculeStack.tif');
end
    

%% display averages image
ns = size(A,3); % selected SM
nr = round(sqrt(ns/2)); % #rows
nc = nr*2; % #columns
A2D = [];
for i = 1:nr
    A2Drow = [];
    for j = 1:nc
        ix = j+(i-1)*nc;
        if ix > ns
            A_ = zeros(wsz);
        else
            A_ = A(:,:,ix);
        end 
        apad = padarray(A_,[1 1]);
        A2Drow = [A2Drow apad];
    end
    A2D = [A2D;A2Drow];
end
imagesc(A2D);
imwrite(uint16(A2D),'singleMoleculeCrops.tif')


colormap('gray')


%% background ratios
reload = 0;
if exist('singleMolData.mat')
    load('singleMolData','RECT')
    reload = 1;
end
ixEx = [4 5 12]; % filter out SM
%ixEx = [1 9 13 14 15 18 20 22];
ixSel(ixEx) = []; 
A(:,:,ixEx) = [];     
ns2 = numel(ixSel);
for i = 1:ns2
    a = A(:,:,i);
    if reload
        [bckgrnd, rect] = imcrop(a,RECT(i,:));
    else
        imagesc(a);
        [bckgrnd, rect] = imcrop;
    end
    peak = a(s-1:s+2,s-1:s+2);
    BCK(i) = mean(bckgrnd(:));
    INT(i) = max(peak(:));
    RECT(i,:) = rect;
end
s = mean(INT);
b = mean(BCK);

(s-b)/b


save('singleMolData','BCK','INT','RECT','ns','ns2','ixEx','ixSel')
    
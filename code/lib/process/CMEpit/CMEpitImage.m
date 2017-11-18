% two parta:
% (1) using selected pits (CMEpitFinder.m) generates a sum recruitment
% image
% (2) using pit crops from time lapse generates a mean image

% E:\MATLAB\SUM-AP2\170521_SUM-AP2\pitselection
clear all
close all

isTL = 0; % 1: time lapse movie; 0: recruitment images 

load('pitCoors.mat'); % LtPits XY
R = imread('binImgRcrtSum.tif');

%% output files
%1 : recs
recPitFN = 'pitRec.tif';
recPitSumFN = 'pitRecSum.tif';
recPitMAT = 'pitRec.mat';
% 2 : TL
selPitFN = 'pitTL.tif';
selPitMeanFN = 'pitTLmean.tif';

ps = 5; % pad size
if ~isTL
    a = 7; % w = 2*a+1
    p = 0;
    for i = 1:size(XY,1)
        Lt = uint16(LtPits(:,:,i));
        y = round(XY(i,2)-0.5)+ps;
        x = round(XY(i,1)-0.5)+ps;
        r = Lt.*R; 
        r = padarray(r,[ps ps]);
        rc = r(y-a+1:y+a,x-a+1:x+a);
        RC(:,:,i) = rc;
        imwrite(rc,recPitFN,'WriteMode','append','Compression', 'none') 
        p = p + rc; % rec sum image
    end
    pitRec = RC;
    save(recPitMAT,'pitRec');
    imwrite(p,recPitSumFN)
    imagesc(p)
    colorbar
end

%% 2 mean TL pit intensity
% run in a folder with processes cropped pits (*_.tif)
% finds the peak intensity frame 
% centers the pits at COM localization

if isTL
    pitsFN = rdir('*_.tif');
    np = numel(pitsFN);
    for i = 1:np
        fn = pitsFN(i).name;
        imginf = imfinfo(fn);
        nf = numel(imginf);
        A=[];
        for j = 1:nf
            A(:,:,j)=imread(fn,j);
        end

        % find highest intensity frame
        Acv = convn(A,ones(3),'same');
        [iy,ix,iz]=ind2sub(size(Acv),find(Acv==max(Acv(:))));
        a = A(iy-4:iy+4,ix-4:ix+4,iz);
        a = repelem(a,4, 4);
        [Xc,Yc] = centOfMassLocCore(a,36);
        dx = round(Xc-18.5);
        dy = round(Yc-18.5);
        a(3+dy:end-2+dy,3+dx:end-2+dx);

        pimg(:,:,i)=a(3+dy:end-2+dy,3+dx:end-2+dx);
    end
    stackWrite(pimg,selPitFN)
    pimgMean = mean(pimg,3);
    imwrite(uint16(pimgMean),selPitMeanFN)
end
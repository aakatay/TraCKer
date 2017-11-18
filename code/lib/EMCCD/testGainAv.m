% using averaged spots
clear; close all;
isCalcResidue = 1;

F = findall(0,'type','figure','tag','TMWWaitbar'); delete(F);
cd C:\MATLAB\TraCKer\data_eclipse\EMCCD

load Spots2 % load spot images
S = S2;
% S = [x,y,frame,spot#]

nSpot = size(S,4);
nAcq = 19;
sp1=1; sp2=nSpot;
sizex = size(S,1);
sizey = sizex;
[X,Y] = meshgrid(double(1:sizex),double(1:sizey));
h2 = waitbar(0,'1- processing for each gain setting...'); 
for i = 1:nAcq
    h = waitbar(0,'2- processing for each spot...'); 
    for j = sp1:sp2
        IMGin = S(:,:,i,j);
        [IMG,DATA] = Gaussian_TIFF(IMGin);
        %[cx,cy,sx,sy,PeakOD,bkgd,pOD] = Gaussian2DJILA(S(:,:,1:10,1))           
        imgDiff = IMG(:,:)-DATA(6) - ( DATA(5)*(exp(-0.5*(((X-DATA(1))/DATA(3)).^2+((Y-DATA(2))/DATA(4)).^2))) );
        res(i,j) = sum(imgDiff(:));
        resMean(i,j) = mean(res(:,i,j),1);
%            sR(:,i) = sqrt(sx(:,i,j).^2+sy(:,i,j).^2); % std
%            ixIgnore = find(sR(:,i) > size(IMG,1) ); % ignore the data
%            sR(ixIgnore,i) = sR(1,i);
%            sRmean(i,j) = mean(sR(:,i));
        waitbar(j / nSpot)
    end
    close(h)
    waitbar(i / nAcq)
end
close(h2)
save('residue','res')
isSaveEachSpot = 1;

if isSaveEachSpot
    for i = 1:nSpot
        %stackWrite(S(:,:,:,i),sprintf('spot%i.tif',i));
    end
end

if isCalcResidue
    for i = sp1:sp2
        figure(2)
        subplot(8,1,i);
        plot(resMean(:,i));
        set(gca,'Xtick',[1:19])
        set(gca,'XTickLabel',['200';'300';'400';'500';'600';'700';'800';'900';'999';'999';'900';'800';'700';'600';'500';'400';'300';'200';'100'])
        xlabel('gain'); ylabel('residue')
    end
else
    for i = sp1:sp2
        figure(2)
        subplot(8,1,i);
        plot(sRmean(:,i));
        set(gca,'Xtick',[1:19])
        set(gca,'XTickLabel',['200';'300';'400';'500';'600';'700';'800';'900';'999';'999';'900';'800';'700';'600';'500';'400';'300';'200';'100'])
        xlabel('gain'); ylabel('std')
    end
end

figure
cc=lines; 
if isCalcResidue
    for i = sp1:sp2
        plot(resMean(:,i),'color',cc(i,:));
        hold on;
    end
    xlabel('gain'); ylabel('residue')
    title('spots: residue vs. gain')
else
    for i = sp1:sp2
        plot(sRmean(:,i),'color',cc(i,:));
        hold on;
    end    
    xlabel('gain'); ylabel('std')
    title('spots: std vs. gain')
end
hold off;
set(gca,'Xtick',[1:19])
set(gca,'XTickLabel',['200';'300';'400';'500';'600';'700';'800';'900';'999';'999';'900';'800';'700';'600';'500';'400';'300';'200';'100'])
legend('spot1','spot2','spot3','spot4','spot5','spot6','spot7','spot8');grid;


figure
selSpots = [2,3,5];
selSpots = 1:8;
if isCalcResidue
    plot(mean(resMean(:,selSpots),2));
    xlabel('gain'); ylabel('residue')
    title('mean of spots: residue vs. gain')
else
    plot(mean(sRmean(:,selSpots),2));
    xlabel('gain'); ylabel('std')
    title('mean of spots: std vs. gain')
end
set(gca,'Xtick',[1:19])
set(gca,'XTickLabel',['200';'300';'400';'500';'600';'700';'800';'900';'999';'999';'900';'800';'700';'600';'500';'400';'300';'200';'100'])

grid


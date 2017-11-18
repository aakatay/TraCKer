% draws the profile across the structure
% SelectedCells\U373 sample struct\150421\cell10\structProf
close all;
clear all;

R=imread('binImgRcrtSum_time0-75s_185X63Y15x16.tif');
A=imread('AVG_lap_185X63Y15x16.tif');

[n,m] = size(R);
mag = 10;
n=n*mag;
m=m*mag;

pos1 = [ 400, 100]; tit1='rec';
pos2 = [ 200, 100]; tit2='int';
figImg = figure('DoubleBuffer','on','Menubar','none','Name',tit1,'NumberTitle','off','Colormap',hot(256),'Position',[pos1 m n]);
axeImg = axes('Parent',figImg,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren','XLim',0.5+[0 m],'YLim',0.5+[0 n]);

% figImg2 = figure('DoubleBuffer','on','Menubar','none','Name',tit1,'NumberTitle','off','Colormap',gray(256),'Position',[pos2 m n]);
% axeImg2 = axes('Parent',figImg2,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren','XLim',0.5+[0 m],'YLim',0.5+[0 n]);

imagesc(R)
axis image
if exist('profStructData.mat')
    load('profStructData');
else
    P1 = [4 8; 8 8];
end
him = imline(gca,[P1(1),P1(2)],[P1(3),P1(4)]);
hold on ; hs=scatter(P1(1),P1(3),'r','o','filled'); 
hold off;
%set(gca,'Xlim',xl)
%set(gca,'Ylim',yl)

P0=getPosition(him);
while 1 % adjust line position
    figure(1)
    waitforbuttonpress
 
    % detect key
    P1=getPosition(him); 
    Pc = sum((P0-P1)~=0,2); % detect change
    Pce = find(Pc>0); % index of the end changed position
    if isempty(Pce) % no changes
        kp=get(gcf,'CurrentCharacter')+1; 
        if kp == 14 % enter key (finish)
            break;
        end 
    end
    
    delete(hs)
    hold on ; hs=scatter(P1(1),P1(3),'r','o','filled'); hold off;
    pause(0.1)
    figStructProfDraw; % e1x e2x,e1y e2y --> prof0, prof2

    P0=getPosition(him);
end

save('profStructData','P1','profInt','profRec'); % 

% print rec image
delete(hs);delete(him)
line([P1(1),P1(2)],[P1(3),P1(4)])
imgZFout = 'profStructLineDraw.tif';
imgFig = getframe(figImg);
imgOut = imgFig.cdata;
imwrite(flipud(imgOut),imgZFout,'Compression', 'none') 

%% print profile plots
figure(101);
imgZFout = 'profStructLinePlots2.tif';
imgFig = getframe;
imgOut = imgFig.cdata;
imwrite(imgOut,imgZFout,'Compression', 'none') 

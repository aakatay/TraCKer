clear
close all;
%clf(201); clf(202); clf(203); clf(204); clf(205); clf(206)
% figure(99)
% y=randn(30,80); x=1:size(y,2);
% shadedErrorBar(x,mean(y,1),std(y),'g');
% shadedErrorBar(x,y,{@median,@std},{'r-o','markerfacecolor','r'});    
% shadedErrorBar([],y,{@median,@std},{'r-o','markerfacecolor','r'}); 

%% print figures

isNormEachProf = 0; % BA ( before or after)
isNormEachBeforeProf = 1; % B
isNormEachBeforeAfterProf = 0; % B-A ( before and after)

alpha = 0.15;
%alpha=.1;
profFN  = 'profDataSel.mat';
load(profFN);

% profile select        
%plen = sum(profC0>=0,2); % profile length
%psel = find(plen>0);
psel = find(sum(exy>0,2)>0); % defined profiles

pxSz=0.16;

D = 11;
xsz = size(profC0,2);
x1 = xsz/2-D/2+1;
x2 = xsz/2+D/2;

np = numel(psel);
exy = exy(psel,:);
pC0 = profC0(psel,x1:x2);
pC2 = profC2(psel,x1:x2);
pCl = profClate(psel,x1:x2);
pD0 = profD0(psel,x1:x2);
pD2 = profD2(psel,x1:x2);   
pDl = profDlate(psel,x1:x2);
plen = sum(profC0(psel,:)>=0,2); % profile length

%% normalize wrt cell
if numel(nProf)>1
    nlast = 0;
    for i = 1:numel(nProf) % cell int coeff
        cellIx = nlast+1:nlast+nProf(i);
        minC = min(min([pC0(cellIx,:); pC2(cellIx,:)]));
        minD = min(min([pD0(cellIx,:); pD2(cellIx,:)]));
        maxC = max(max([pC0(cellIx,:); pC2(cellIx,:)]));
        maxD = max(max([pD0(cellIx,:); pD2(cellIx,:)]));
        cellIntensity(i) = maxC-minC;
    %         pC0(cellIx,:) = (pC0(cellIx,:)-minC)/(maxC-minC);
    %         pD0(cellIx,:) = (pD0(cellIx,:)-minD)/(maxD-minD);
    %         pC2(cellIx,:) = (pC2(cellIx,:)-minC)/(maxC-minC);
    %         pD2(cellIx,:) = (pD2(cellIx,:)-minD)/(maxD-minD);
        nlast = cellIx(end);
    end

    refCellIx = 2;
    cellIntensityCoeff = cellIntensity/cellIntensity(2);

    nlast = 0;
    for i = 1:numel(nProf)
        cellIx = nlast+1:nlast+nProf(i);
        pC0(cellIx,:) = minC+(pC0(cellIx,:)-minC)/cellIntensityCoeff(i);
        pD0(cellIx,:) = minD+(pD0(cellIx,:)-minD)/cellIntensityCoeff(i);
        pC2(cellIx,:) = minC+(pC2(cellIx,:)-minC)/cellIntensityCoeff(i);
        pD2(cellIx,:) = minD+(pD2(cellIx,:)-minD)/cellIntensityCoeff(i);
        nlast = cellIx(end);
    end
    % psel=142;
    % exy(psel,:)=[];
    % profC0(psel,:)=[];
    % profC2(psel,:)=[];
    % profClate(psel,:)=[];
    % profD0(psel,:)=[];
    % profD2(psel,:)=[];
    % profDlate(psel,:)=[];
end

%% interp coors
dp = linspace(-D/2,D/2,D*10+1);
for i = 1:np
    pl = plen(i);
    px=(pl-D)/2;
    pC0(i,:);
    x1 = exy(i,1);
    x2 = exy(i,3);
    xm = (x1+x2)/2;
    x1 = x1-xm; x2 = x2-xm;
    x=linspace(x1,x2,pl);
    x=x(px+1:end-px);
    xl=x(end)-x(1);
    y1 = exy(i,2);
    y2 = exy(i,4);
    ym = (y1+y2)/2;
    y1 = y1-ym; y2 = y2-ym;
    y=linspace(y1,y2,pl);
    y=y(px+1:end-px);
    yl=y(end)-y(1);
    d = sqrt(x.^2+y.^2);
    d(1:5)=-d(1:5);
    pC0p(i,:) = interp1(d,pC0(i,:),dp);
    pC2p(i,:) = interp1(d,pC2(i,:),dp);
    pClp(i,:) = interp1(d,pCl(i,:),dp);
    pD0p(i,:) = interp1(d,pD0(i,:),dp);
    pD2p(i,:) = interp1(d,pD2(i,:),dp);
    pDlp(i,:) = interp1(d,pDl(i,:),dp);
    %figure(89); plot(pC0p(i,:)); figure(839); plot(pC0(i,:))
    ccc=4;
end
lenp = sum(~isnan(pC0p),2); %lengths in interpolated coors
ixnan=(D*10+1-min(lenp))/2;
pC0p = pC0p(:,ixnan+1:end-ixnan);
pC2p = pC2p(:,ixnan+1:end-ixnan);
pClp = pClp(:,ixnan+1:end-ixnan);
pD0p = pD0p(:,ixnan+1:end-ixnan);
pD2p = pD2p(:,ixnan+1:end-ixnan);
pDlp = pDlp(:,ixnan+1:end-ixnan);

dp = dp(ixnan+1:end-ixnan);
dl = numel(dp);
vx = pxSz*dp;

%% normalization

% scalar normalization
minC = min( [pC0p(:); pC2p(:)] );
maxC = max( [pC0p(:); pC2p(:)] );

minD = min( [pD0p(:); pD2p(:)] );
maxD = max( [pD0p(:); pD2p(:)] );


% vector normalization (wrt each profile)
minC0 = min(pC0p,[],2); 
minC2 = min(pC2p,[],2); 
maxC0 = max(pC0p,[],2); 
maxC2 = max(pC2p,[],2); 
minD0 = min(pD0p,[],2); 
minD2 = min(pD2p,[],2); 
maxD0 = max(pD0p,[],2); 
maxD2 = max(pD2p,[],2); 
if isNormEachProf
    minC = min([minC0 minC2],[],2);
    maxC = max([maxC0 maxC2],[],2);
    minD = min([minD0 minD2],[],2);
    maxD = max([maxD0 maxD2],[],2);
    minC = [minC; minC];
    maxC = [maxC; maxC];
    minD = [minD; minD];
    maxD = [maxD; maxD];
elseif isNormEachBeforeProf
    minC = [minC0; minC0];
    maxC = [maxC0; maxC0];
    minD = [minD0; minD0];
    maxD = [maxD0; maxD0];
elseif isNormEachBeforeAfterProf
    minC = [minC0; minC2];
    maxC = [maxC0; maxC2];
    minD = [minD0; minD2];
    maxD = [maxD0; maxD2];
end
if isNormEachProf | isNormEachBeforeProf | isNormEachBeforeAfterProf
    minC = repmat(minC,1,dl);
    maxC = repmat(maxC,1,dl);
    minD = repmat(minD,1,dl);
    maxD = repmat(maxD,1,dl);
end

% array
pcComp = [pC0p; pC2p];
pdComp = [pD0p; pD2p];

Dx = 81;
xLim=[1:Dx]+ (dl-Dx)/2;
pcComp0 = pcComp(1:np,xLim);
pcComp2 = pcComp(1+np:end,xLim);
pdComp0 = pdComp(1:np,xLim);
pdComp2 = pdComp(1+np:end,xLim);
pcComp02 = pcComp0-pcComp2;
pdComp02 = pdComp0-pdComp2;
mpcComp = [mean(pcComp0,1); mean(pcComp2,1); mean(pcComp02,1)]; % means
mpdComp = [mean(pdComp0,1); mean(pdComp2,1); mean(pdComp02,1)];

% normalized
pcCompNorm = (pcComp-minC)./(maxC-minC);
pdCompNorm = (pdComp-minD)./(maxD-minD);

pcCompNorm0 = pcCompNorm(1:np,xLim);
pcCompNorm2 = pcCompNorm(1+np:end,xLim);
pdCompNorm0 = pdCompNorm(1:np,xLim);
pdCompNorm2 = pdCompNorm(1+np:end,xLim);
vx = vx(xLim);


pcCompNorm02 = pcCompNorm0-pcCompNorm2; % mpC0-mpC2
pdCompNorm02 = pdCompNorm0-pdCompNorm2; % mpD0-mpD2

% normalized mean
pcCompNormMean = [mean(pcCompNorm0,1); mean(pcCompNorm2,1); mean(pcCompNorm02,1)]'; % means
pdCompNormMean = [mean(pdCompNorm0,1); mean(pdCompNorm2,1); mean(pdCompNorm02,1)]';

%%
C1 = [0 0 1];
C2 = [1 0.5 0];
C3 = [0 0.8 0];
lt0=6; lt=3;

%%
figure(201); 
plot(vx-min(vx),mpcComp); %title(sprintf('clathrin mean N=%i',np)); 
%hl=legend('C0','C1','C0-C1');set(hl,'EdgeColor',[1 1 1])
xlabel('distance [\mum]')
ylabel('intensity (au)')
axis tight
yticks=get(gca,'YTick');
yticks_ = [0 yticks(end) ];
set(gca,'YTick',yticks_); 
set(gcf,'color','w');
set(gca,'FontSize',15)
gcaChild = get(gca,'Children');
set(gcaChild(3),'Color',C1,'LineWidth',lt0)
set(gcaChild(2),'Color',C2,'LineWidth',lt0)
set(gcaChild(1),'Color',C3,'LineWidth',lt0)
linespec= {'Color' [1 0.5 0]};
set(gca, 'box', 'off')

%%
figure(202);  clf(202)
hs=shadedErrorBar(vx,pcCompNormMean(:,1),std(pcCompNorm0),'b',alpha);
hold on
hs=shadedErrorBar(vx,pcCompNormMean(:,2),std(pcCompNorm2),linespec,alpha);
%plot(vx,pcCompNorm(:,2));
hold off
%hl=legend('C0','C1','C0-C1');set(hl,'EdgeColor',[1 1 1])
xlabel('distance [\mum]')
%ylabel('intensity (au)')
%xticks=get(gca,'XTick');
xticks = [-0.64 0.64];
xticks_ = [xticks(1) xticks(1)/2 0 xticks(end)/2 xticks(end)];
set(gca,'XTick',xticks_);
yticks=get(gca,'YTick');
yticks_ = [0:0.2:yticks(end)];
set(gca,'YTick',yticks_);
set(gcf,'color','w');
set(gca,'FontSize',15)
gcaChild = get(gca,'Children');
%set(gcaChild(2),'Color',C1,'LineWidth',lt)
%set(gcaChild(1),'Color',C2,'LineWidth',lt)
set(gcaChild(5),'LineWidth',lt)
set(gcaChild(1),'LineWidth',lt)
axis tight
yLim=ylim;
set(gca,'YLim',[yLim(1) yticks_(end)])

set(gca, 'box', 'off','xcolor','w')

%%
figure(203); clf(203)
shadedErrorBar(vx,pcCompNormMean(:,3),std(pcCompNorm02),'g'); %title(sprintf('clathrin mean N=%i',np)); 
%hl=legend('C0','C1','C0-C1');set(hl,'EdgeColor',[1 1 1])
xlabel('distance [\mum]')
%ylabel('intensity (au)')
%xticks=get(gca,'XTick');
xticks = [-0.64 0.64];
xticks_ = [xticks(1) xticks(1)/2 0 xticks(end)/2 xticks(end)];
set(gca,'XTick',xticks_);
yticks=get(gca,'YTick');
yticks_ = [yticks(1):0.2:yticks(end)];
set(gcf,'color','w');
set(gca,'FontSize',15)
gcaChild = get(gca,'Children');
set(gcaChild(1),'LineWidth',lt)
axis tight
yLim=ylim;
set(gca,'YLim',[yLim(1) yticks_(end)])

set(gca,'YTick',yticks_);
set(gca, 'box', 'off')
set(gca,'YLim',[yLim(1) yticks_(end)])


%% ================== dynamin ==================

figure(204); 
plot(vx-min(vx),mpdComp); %title(sprintf('clathrin mean N=%i',np)); 
%hl=legend('C0','C1','C0-C1');set(hl,'EdgeColor',[1 1 1])
xlabel('distance [\mum]')
ylabel('intensity (au)')
axis tight
yticks=get(gca,'YTick');
yticks_ = [0 yticks(end) ];
set(gca,'YTick',yticks_);
set(gcf,'color','w');
set(gca,'FontSize',15)
gcaChild = get(gca,'Children');
set(gcaChild(3),'Color',C1,'LineWidth',lt0)
set(gcaChild(2),'Color',C2,'LineWidth',lt0)
set(gcaChild(1),'Color',C3,'LineWidth',lt0)
set(gca, 'box', 'off')
%%
figure(205); clf(205)
hs=shadedErrorBar(vx,pdCompNormMean(:,1),std(pcCompNorm0),'b',alpha);
hold on
hs=shadedErrorBar(vx,pdCompNormMean(:,2),std(pcCompNorm2),linespec,alpha);
%plot(vx,pdCompNorm(:,2)); %title(sprintf('clathrin mean N=%i',np)); 
%hl=legend('C0','C1','C0-C1');set(hl,'EdgeColor',[1 1 1])
xlabel('distance [\mum]')
%ylabel('intensity (au)')
%xticks=get(gca,'XTick');
xticks = [-0.64 0.64];
xticks_ = [xticks(1) xticks(1)/2 0 xticks(end)/2 xticks(end)];
set(gca,'XTick',xticks_);set(gca, 'box', 'off')
yticks=get(gca,'YTick');
yticks_ = [0:0.2:yticks(end)];
set(gca,'YTick',yticks_);
set(gcf,'color','w');
set(gca,'FontSize',15)
gcaChild = get(gca,'Children');
%set(gcaChild(2),'Color',C1,'LineWidth',lt)
%set(gcaChild(1),'Color',C2,'LineWidth',lt)
set(gcaChild(5),'LineWidth',lt)
set(gcaChild(1),'LineWidth',lt)
set(gca,'XTick',[]);
axis tight
set(gca, 'box', 'off','xcolor','w')

%%
figure(206); 
shadedErrorBar(vx,pdCompNormMean(:,3),std(pcCompNorm02),'g');%title(sprintf('clathrin mean N=%i',np)); 
%hl=legend('C0','C1','C0-C1');set(hl,'EdgeColor',[1 1 1])
xlabel('distance [\mum]')
%ylabel('intensity (au)')
%xticks=get(gca,'XTick');
xticks = [-0.64 0.64];
xticks_ = [xticks(1) xticks(1)/2 0 xticks(end)/2 xticks(end)];
set(gca,'XTick',xticks_);
yticks=get(gca,'YTick');
yticks_ = [0:0.2:yticks(end)];
set(gca,'YTick',yticks_);
set(gcf,'color','w');
set(gca,'FontSize',15)
gcaChild = get(gca,'Children');
set(gcaChild(1),'LineWidth',lt)
axis tight
set(gca, 'box', 'off')

%% adjust scales

figure(202); 
ylim = get(gca,'Ylim');
posBA = get(gca,'Position'); % before and after
ylBA = ylim(2)-ylim(1);
figure(203); 
ylim = get(gca,'Ylim');
posDiff = get(gca,'Position'); % diff
ylDiff = ylim(2)-ylim(1);

posDiff(4)=posBA(4)/ylBA*ylDiff;
set(gca,'Position',posDiff);

figure(205); 
ylim = get(gca,'Ylim');
posBA = get(gca,'Position'); % before and after
ylBA = ylim(2)-ylim(1);
figure(206); 
ylim = get(gca,'Ylim');
posDiff = get(gca,'Position'); % diff
ylDiff = ylim(2)-ylim(1);

posDiff(4)=posBA(4)/ylBA*ylDiff;
set(gca,'Position',posDiff);


%%

figure(201); 
xlim0 = get(gca,'Xlim');
xlim = xlim0+[-0.05 0.05];
set(gca,'xlim',xlim);

set(gca,'linewidth',1.5)

figure(202); 
set(gca,'linewidth',2)
figure(203); 
set(gca,'linewidth',2)




%% print figures
figno=[201:206];
imgZFout = 'figs.tif';
delete(imgZFout);
for i = 1:numel(figno)
    figure(figno(i));
    
    imgFig = getframe(gcf); 
    imgOut = imgFig.cdata;
    imwrite(imgOut,imgZFout,'WriteMode','append','Compression', 'none') 

end
close(204)
close(205)
close(206)

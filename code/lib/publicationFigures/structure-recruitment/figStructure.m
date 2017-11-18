% run findStructs.m before and select structure of interest then run this


isSaveFigs = 1;

if ~isSaveFigs    % load
    ixS=38;
    cellLabel = '160726-cell7';
    ixS=12;
    cellLabel = '150421-cell10';
else
    dirName = pwd; % directory name
    dS = strfind(dirName,'\'); % dirSlash
    cellName = dirName(dS(end-1)+1:dS(end)-1);
    dateName = dirName(dS(end-2)+1:dS(end-2)+6);
    cellLabel = [dateName '-' cellName maskLab];
end


fig1FN = sprintf('%s-str%03i_fig1.fig',cellLabel,ixS);
fig2FN = sprintf('%s-str%03i_fig2.fig',cellLabel,ixS);
fig3FN = sprintf('%s-str%03i_fig3.fig',cellLabel,ixS);
fig4FN = sprintf('%s-str%03i_fig4.fig',cellLabel,ixS);
fig10FN = sprintf('%s-str%03i_fig10.fig',cellLabel,ixS);
fig999FN = sprintf('%s-str%03i_fig999.fig',cellLabel,ixS);

if isSaveFigs
    figure(2)
    savefig(fig2FN)
    figure(3)
    savefig(fig3FN)
    figure(4)
    savefig(fig4FN)
    figure(10)
    savefig(fig10FN)
    figure(999)
    savefig(fig999FN)
end
    
close all;
openfig(fig2FN); 
openfig(fig3FN); 
openfig(fig4FN);
openfig(fig10FN);
openfig(fig999FN);

imgZFout = sprintf('%s-%03i_struct.tif',cellLabel,ixS);
imgZFout2 = sprintf('%s-%03i_structProf.tif',cellLabel,ixS);
imgZFout3 = sprintf('%s-%03i_colorBars.tif',cellLabel,ixS);
delete(imgZFout);
delete(imgZFout3);


%% process border extension
figure(1);% recs
GCA=gca;
GCALine=GCA.Children;
ds=0;
set(GCALine(1+ds),'Visible','off')
R = GCALine(2+ds).CData;
xl = GCALine(1+ds).XData;
yl = GCALine(1+ds).YData;

b=[yl' xl'];
szXY = size(R);
[b]=extendBoundaries(b,szXY);


nr = max(R(:));
CMfull = hot(64);
cix = round(linspace(1,64,nr+1));
CMdiscrete = CMfull(cix,:);
colormap(CMdiscrete);

hold on;
plot(b(:,2),b(:,1),'g','LineWidth',1.5)
hold off;


%% save (1)
figure(1)

imgFig = getframe(gcf); 
imgOut = imgFig.cdata;
imwrite(imgOut,imgZFout,'Compression', 'none') 

%% process
figure(2);% profiles

figure(3); % prebleach
GCA=gca;
GCALine=GCA.Children;
set(GCALine(1),'Visible','off')
A = GCALine(2).CData;
figure(5); % mask
GCA=gca;
GCALine=GCA.Children;
%set(GCALine(1),'Visible','off')
M = GCALine(1).CData;

%% save (2)
figure(5); 
CM1=[0 0 0; 1 0.35 0]; colormap(CM1);
imgFig = getframe(gcf); 
imgOut = imgFig.cdata;
imwrite(imgOut,imgZFout,'WriteMode','append','Compression', 'none') 

%% save(3)
figure(3); % prebleach
load mx
A=A/64*mx;
A = M.*A;
A(A==0)=min(A(A~=0))-max(A(:)/15);
imagesc(A);
CM=colormap('hot');
CM(end-8:end,:)=[];
colormap(CM)
imgFig = getframe(gcf); 
imgOut = imgFig.cdata;
imwrite(imgOut,imgZFout,'WriteMode','append','Compression', 'none') 

%% save(prof)
figure(2); 
posCtrl = [100 200 800 600];
set(gcf,'units','pixels','Position',posCtrl);
set(gca,'units','pixels','Position',posCtrl+[0 -100 -160 -170]); 

title('');
GCA=gca;
Lines=GCA.Children;
th = 4;
set(Lines(1),'LineWidth',th,'Color',[1 0 0])
set(Lines(2),'LineWidth',th,'Color',[0 0.5 0])
set(Lines(3),'LineWidth',th,'Color',[0 0 1])

px = 160/1.5/4;
xx=[15:50:400]; % [nm]
xd = (xx - px/2)/px;

xticks=get(gca,'XTick');
nt = sum(xd<xticks(end)); % number of ticks
xticks_ = xd(1:nt);
set(gca,'XTick',xticks_);

yticks=get(gca,'YTick');
yt = round(yticks(end)/5/20)*20;
yticks_ = [0:yt:yticks(end)];
yticks_(end+1) = yticks_(end)+yt;
set(gca,'YTick',yticks_);

set(gca,'linewidth',2)
set(gca,'TickDir','out')
grid off
set(gca, 'box', 'off')

xlim0 = get(gca,'Xlim');
xlim = xlim0+[-0.5 0.5];
set(gca,'xlim',xlim);

ylim0 = get(gca,'Ylim');
ylim = ylim0+[-10 10];
set(gca,'ylim',ylim);

set(gcf,'color','w');
hb=findall(0,'Type','UIControl'); delete(hb);
imgFig = getframe(gcf); 
imgOut = imgFig.cdata;
imwrite(imgOut,imgZFout2,'Compression', 'none') 


%legend('recruitment','uniform','intensity')

%% colorbars
delete(imgZFout3);
figure(3)
nc = int2str([0:nr]');
for i=[1 3 5]
    figure(i);
    
    %figure; imagesc(R)  

     b = get(gca,'children');% link to the data curve/s
     % image data
    figure(777)
    
    A = b(end).CData;
    
    if i == 1
        colormap(CMdiscrete);
    elseif i == 3 % intensity image
        colormap(CM);
    elseif i == 5 % intensity image
        colormap(CM1);
    end
    imagesc(A)
    hc=colorbar;
    if i == 1, set(hc,'YTick',linspace(0.4,nr-0.4,nr+1),'YTickLabel',nc); end


    
     a = get(gcf,'children');% link to ledgend is now a(1)
     b = get(gca,'children');% link to the data curve/s
     
     set(a(2),'visible','off'); %hide axes etc...
%      set(b,'visible','off'); %hide data...
%      set(b,'visible','on'); %hide data...

     legfs = get(a(1),'Fontsize'); %get legend fontsize
     set(a(1),'Fontsize',legfs+10); %make legend appear larger

    set(gcf,'color','w');
    set(gca,'units','pixels'); 
    p = get(gca,'Position');
    p(4)=700/2;
    set(gca,'units','pixels','Position',[-p(3) 25 p(3) p(4)]); 

    p = get(gcf,'Position');
    p(4)=850/2;
    set(gcf,'units','pixels','Position',[100 100 100 p(4)]); 
    hb=findall(0,'Type','UIControl'); delete(hb);
    
    imgFig = getframe(gcf); 
    imgOut = imgFig.cdata;
    imwrite(imgOut,imgZFout3,'WriteMode','append','Compression', 'none') 
    delete(777)

end
%close all











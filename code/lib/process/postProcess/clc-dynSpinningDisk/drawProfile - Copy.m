% draws profile for profile selection (in selectStruct)
figure(97); % subplot images
gsb=(findall(gcf,'type','axes')); % sub plots

nsb=length(gsb); % number of subplots
isDel=0;
if exist('hl'), isDel=1; end;
for isb = 1:nsb
    axes(gsb(isb))
    hold on;
    if isDel, delete(hl(isb)); end;
    hl(isb)=line([e1x e2x],[e1y e2y],'Color','r');
    xAll=[]; yAll=[];
    if strcmp(get(get(gca,'title'),'String'),'Adiff(t+1)-(t)')
        xAll = xc1All;
        yAll = yc1All;
    elseif strcmp(get(get(gca,'title'),'String'),'Adiff2(t+2)-(t+1)')
        xAll = xc2All;
        yAll = yc2All;
    elseif strcmp(get(get(gca,'title'),'String'),'Ddiff(t+1)-(t)')
        xAll = xd1All;
        yAll = yd1All;
    elseif strcmp(get(get(gca,'title'),'String'),'Ddiff2(t+2)-(t+1)')
        xAll = xd2All;
        yAll = yd2All;
    end    
	scatter(xAll,yAll,'k','.')
    
    scatter(xd0,yd0,'r','.')
    scatter(xc0,yc0,'c','.')
    hold off;
    set(gca,'Xlim',xl)
    set(gca,'Ylim',yl)
end


%% profiles

[px1,py1,p0_1] = improfile(A0_2,[e1x xd0],[e1y yd0]);
[px2,py2,p0_2] = improfile(A0_2,[xd0 e2x],[yd0 e2y]);
p0 = [p0_1(1:end-1); p0_2]; % merge
px = [px1(1:end-1); px2]; % data point coors
py = [py1(1:end-1); py2];
ixC = max(find(p0==A0_2(round(yd0),round(xd0)))); %center index
p0(isnan(p0)) = [];
pl = numel(p0);
if ~rem(pl,2)
    warning('not centered, skipping')
    %continue
end
    
p2 = improfile(A2,px,py);
p3 = improfile(A3,px,py);
p4 = improfile(A4,px,py);
pLate = improfile(Alate,px,py);
p0d = improfile(D0_2,px,py);
p2d = improfile(D2,px,py);
pdLate = improfile(Dlate,px,py);

p2(isnan(p2)) = [];
p3(isnan(p3)) = [];
p4(isnan(p4)) = [];
pLate(isnan(pLate)) = [];
p0d(isnan(p0d)) = [];
p2d(isnan(p2d)) = [];
pc = floor((pl0-pl)/2)+1;
pdLate(isnan(pdLate)) = [];

if ixC~=ceil(pl/2), 
    %pc = pc -1;
end

p2c = p2;
if tc~=td
    p2c = p3;
end

profC0(i,:) = nan;
profC2(i,:) = nan;
profClate(i,:) = nan;
profD0(i,:) = nan;
profD2(i,:) = nan;
profDlate(i,:) = nan;
    
profC0(i,pc:pc+pl-1) = p0;
profC2(i,pc:pc+pl-1) = p3;
profClate(i,pc:pc+pl-1) = pLate;
profD0(i,pc:pc+pl-1) = p0d;
profD2(i,pc:pc+pl-1) = p2d;
profDlate(i,pc:pc+pl-1) = pdLate;

if isDisp, return; end;
% clathrin profiles before and after dyn peak
figure(101) % profile
sz0 = 1.5; sz1 = 1.5; sz2 = 0.75; sz3 = 0.75;
if tc~=td, sz2=1.5; sz1=0.75; end;
plot(p0,'k','LineWidth',sz0);
hold on;
plot(p2,'b','LineWidth',sz1); % 1-1 frame after
plot(p3,'g','LineWidth',sz2); % 2-2
plot(p4,'r','LineWidth',sz3); % 3-3
plot(p0d-p2d,'c','LineWidth',sz3); % 3-3
hold off
line([ixC ixC],[0 max([p0' p2'])])
legend('t-1:t', 't+1', 't+2', 't+3','dyn')
title(sprintf('clathrin intensity line profile, peak frames: #%i & #%i (clathrin&dynamin)',tc,td))
grid minor;
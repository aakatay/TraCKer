% search and compile results from different cells
%close all


sz=[200 400]; % [px] small mid large separation

ixDel = [6 11 12 13];
ixDel=[];
matFNdir = rdir('**\edgeDistData*.mat');


save('matFNdir','matFNdir')
edgeDistMeanArrAllCells = [];
ncells = numel(matFNdir);
for i = 1: ncells
    matFN = matFNdir(i).name;
    ixSl = strfind(matFN,'\');
    cellIx = str2num(matFN(5:ixSl-1));
    if ismember(cellIx,ixDel), disp('skipping cell');ncells=ncells-1; continue; end;
    load(matFN,'edgeDistMeanArr');
    edgeDistMeanArrAllCells = [edgeDistMeanArrAllCells; edgeDistMeanArr];
end
edgeDistMeanArr = edgeDistMeanArrAllCells;


%% plot scatter maps : dist vs area
figure(100);
x = edgeDistMeanArr(:,4); % area
y = edgeDistMeanArr(:,1)./edgeDistMeanArr(:,2); % vs uniform recruitment
z = edgeDistMeanArr(:,1)./edgeDistMeanArr(:,3); % vs intensity

scatter(x,y,'.')    
mxx = max(x);
mxy = max(y);
nAbove = sum(y>1);
nBelow = sum(y<=1);
N = sum(y>0);
ratAbove = nAbove/N;
ratBelow = nBelow/N;
line([0 mxx], [1 1]);
text(mxx*0.8,mxy*0.9,sprintf('>1:\t%.02f (%i)\n<=1:\t%.02f (%i)',ratAbove,nAbove,ratBelow,nBelow));
yl=ylim;
line([sz(1) sz(1)],[yl(1) yl(2)],'Color','r');
line([sz(2) sz(2)],[yl(1) yl(2)],'Color','r');

xlabel('area');
ylabel('edge distance (norm.)')
grid minor


fnameXLS=sprintf('edgeDistData-%icells.xls',ncells);

fnameEdgeDist=sprintf('edgeDistvsAreaScatter-normUniform-%icells.tif',ncells);

fnameEdgeDistHist1=sprintf('edgeDistRatAreaBlocks-normUniform-%icells.tif',ncells);
fnameEdgeDistHist2=sprintf('edgeDistRatAreaBlocks-normIntensity-%icells.tif',ncells);



fnameEdgeDistHistvsArea{1}=sprintf('edgeDistRatHist-normUniform_All-%icells.tif',ncells);
fnameEdgeDistHistvsArea{2}=sprintf('edgeDistRatHist-normUniform_Sma-%icells.tif',ncells);
fnameEdgeDistHistvsArea{3}=sprintf('edgeDistRatHist-normUniform_Mid-%icells.tif',ncells);
fnameEdgeDistHistvsArea{4}=sprintf('edgeDistRatHist-normUniform_Lar-%icells.tif',ncells);

fnameEdgeDistHistvsArea{5}=sprintf('edgeDistRatHist-normIntnsty_All-%icells.tif',ncells);


title(sprintf('edge distance -norm. by Uniform- (N = %i)',N));
imgFig = getframe(gcf);
dataImg = imgFig.cdata; 
imwrite(uint16(dataImg),fnameEdgeDist); % structMap

%% stacked area blocks 
figure(101)
y1 = y(x<=sz(1));
x1 = x(x<=sz(1));
y2 = y(find((sz(1)<x).*(x<=sz(2))));
x2 = x(find((sz(1)<x).*(x<=sz(2))));
y3 = y(sz(2)<x);
x3 = x(sz(2)<x);
n1=numel(y1);
n2=numel(y2);
n3=numel(y3);
nmx = max([n1 n2 n3]);
yp = nan(nmx,3);
yp(1:n1,1)=y1;
yp(1:n2,2)=y2;
yp(1:n3,3)=y3;
[br1,bins] = histcounts(yp(:,1)',[0 1 1000]);
[br2,bins] = histcounts(yp(:,2)',[0 1 1000]);
[br3,bins] = histcounts(yp(:,3)',[0 1 1000]);
[brAll,bins] = histcounts(yp(:)',[0 1 1000]);
br = [brAll; br1 ; br2 ; br3];
br = br./repmat(sum(br,2),1,size(br,2))*100;
%bar([0.5 1.5],br');
bar(br,'stacked')
set(gca,'XTick',([1:4]));
set(gca,'XTickLabels',({'all','small','mid','large'}))
ltx1 = sprintf('distRatio<=1');
ltx2 = sprintf('distRatio>1');
legend(gca,ltx1,ltx2,'Location','northwest')
legend(gca,ltx1,ltx2,'Location','northeast')
xlabel('edge distance (norm.)')
ylabel('distribution (percentage)')
title(sprintf('edge distance histogram -norm. by Uniform- (N = %i)',N));
grid minor
imgFig = getframe(gcf);
dataImg = imgFig.cdata; 
imwrite(uint16(dataImg),fnameEdgeDistHist1); % structMap
yBR = [brAll;br1;br2;br3];

figure(201)
z1 = z(x<=sz(1));
z2 = z(find((sz(1)<x).*(x<=sz(2))));
z3 = z(sz(2)<x);
nz1=numel(z1);
nz2=numel(z2);
nz3=numel(z3);
nmx = max([nz1 nz2 nz3]);
zp = nan(nmx,3);
zp(1:nz1,1)=z1;
zp(1:nz2,2)=z2;
zp(1:nz3,3)=z3;
[br1,bins] = histcounts(zp(:,1)',[0 1 1000]);
[br2,bins] = histcounts(zp(:,2)',[0 1 1000]);
[br3,bins] = histcounts(zp(:,3)',[0 1 1000]);
[brAll,bins] = histcounts(zp(:)',[0 1 1000]);
br = [brAll; br1 ; br2 ; br3];
br = br./repmat(sum(br,2),1,size(br,2))*100;
%bar([0.5 1.5],br');
bar(br,'stacked')
set(gca,'XTick',([1:4]));
set(gca,'XTickLabels',({'all','small','mid','large'}))
ltx1 = sprintf('distRatio<=1');
ltx2 = sprintf('distRatio>1');
legend(gca,ltx1,ltx2,'Location','northwest')
legend(gca,ltx1,ltx2,'Location','northeast')
xlabel('edge distance (norm.)')
ylabel('distribution (percentage)')
title(sprintf('edge distance histogram -norm. by Intensity- (N = %i)',N));
grid minor
imgFig = getframe(gcf);
dataImg = imgFig.cdata; 
imwrite(uint16(dataImg),fnameEdgeDistHist2); % structMap.
zBR = [brAll;br1;br2;br3];

%% area histograms
%figure;
mxrat = ceil(max(yp(:))*10)/10; % max ratio
mnrat = floor(min(yp(:))*10)/10; % min ratio
db=0.05;
histbins = mnrat:db:mxrat;
histbins = histbins-db/2;
histbins = histbins(2:end);

mxrat = ceil(max(zp(:))*10)/10; % max ratio
mnrat = floor(min(zp(:))*10)/10; % min ratio
db=0.05;
histbins2 = mnrat:db:mxrat;
histbins2 = histbins2-db/2;
histbins2 = histbins2(2:end);


N = numel(y);
p1 = round(n1/N*100);
p2 = round(n2/N*100);
p3 = round(n3/N*100);
ltx = sprintf('all structures, N=%i',N);
ltx1 = sprintf('area<=%i [px], N=%i(%i%%)',sz(1),n1,p1);
ltx2 = sprintf('%i<area<=%i [px], N=%i(%i%%)',sz(1),sz(2),n2,p2);
ltx3 = sprintf('area>%i [px], N=%i(%i%%)',sz(2),n3,p3);

%% area histograms figures
% 1: all
figure(11)
hist(yp(:),histbins)
line([1 1],[0 max(hist(yp(:),histbins))],'Color','r')
title(ltx)
rate = round(sum(yp(:)<=1)/sum(yp(:)>0)*100);
legend(gca,sprintf('%i%% <= 1',rate),'Location','northeast')
grid minor
xlabel('edge distance (uniform norm.)')
ylabel('# of structures')
% 2: area<=100
figure(12)
hist(yp(:,1),histbins)
line([1 1],[0 max(hist(yp(:,1),histbins))],'Color','r')
title(ltx1)
rate = round(sum(yp(:,1)<=1)/sum(yp(:,1)>0)*100);
legend(gca,sprintf('%i%% <= 1',rate),'Location','northeast')
grid minor
xlabel('edge distance (uniform norm.)')
ylabel('# of structures')
% 3: 100<area<=300
figure(13)
hist(yp(:,2),histbins)
line([1 1],[0 max(hist(yp(:,2),histbins))],'Color','r')
sum(yp(:,1)<=1)/sum(yp(:,1)>0);
title(ltx2)
rate = round(sum(yp(:,2)<=1)/sum(yp(:,2)>0)*100);
legend(gca,sprintf('%i%% <= 1',rate),'Location','northeast')
grid minor
xlabel('edge distance (uniform norm.)')
ylabel('# of structures')

% 4: 300<area
figure(14)
hist(yp(:,3),histbins)
line([1 1],[0 max(hist(yp(:,3),histbins))],'Color','r')
title(ltx3)
rate = round(sum(yp(:,3)<=1)/sum(yp(:,3)>0)*100);
legend(gca,sprintf('%i%% <= 1',rate),'Location','northeast')
grid minor
xlabel('edge distance (uniform norm.)')
ylabel('# of structures')

% 5: all
figure(15)
hist(zp(:),histbins2)
line([1 1],[0 max(hist(zp(:),histbins2))],'Color','r')
title(ltx)
rate = round(sum(zp(:)<=1)/sum(zp(:)>0)*100);
legend(gca,sprintf('%i%% <= 1',rate),'Location','northeast')
grid minor
xlabel('edge distance (intesity norm.)')
ylabel('# of structures')

for i = 1:5
    figure(10+i);
    imgFig = getframe(gcf);
    dataImg = imgFig.cdata; 
    imwrite(uint16(dataImg),fnameEdgeDistHistvsArea{i}); % structMap
end

%% generate XLS cells
xlsShtNM{1} = 'edgeDistRatAreaBlocks-normUniform';
xlsShtNM{2} = 'edgeDistRatAreaBlocks-normIntensity';
xlsShtNM{3} = 'edgeDistRatHist-normUniform_All';
xlsShtNM{4} = 'edgeDistRatHist-normIntnsty_All';


% #1
clls = [{' '} {'all'} {sprintf('area<=%i',sz(1))} {sprintf('%i<area<=%i',sz)} {sprintf('%i<area',sz(2))}]';
clls(2:5,2:3)=mat2cell(yBR,[1 1 1 1],[1 1]);
clls(1,2:3) = [{'<=1'} {'>1'}];
CELLS{1}=clls;

% #2
clls = [{' '} {'all'} {sprintf('area<=%i',sz(1))} {sprintf('%i<area<=%i',sz)} {sprintf('%i<area',sz(2))}]';
clls(2:5,2:3)=mat2cell(zBR,[1 1 1 1],[1 1]);
clls(1,2:3) = [{'<=1'} {'>1'}];
CELLS{2}=clls;

% #3
figure(11)
X_output = get(get(gca,'Children'),'XData');
Y_output = get(get(gca,'Children'),'YData');
X_Data = cell2mat(Y_output(2));
X_Data = X_Data(3,:);


clls = [histbins(1:numel(X_Data)); X_Data];
CELLS{3}=clls;

% #4
figure(15)
X_output = get(get(gca,'Children'),'XData');
Y_output = get(get(gca,'Children'),'YData');
X_Data = cell2mat(Y_output(2));
X_Data = X_Data(3,:);


clls = [histbins2(1:numel(X_Data)); X_Data];
CELLS{4}=clls;

%% write to XLS
delete(fnameXLS);

% index
shtNM = [{'2'} {'3'} {'4'} {'5'}]';
shtNM(1:4,2) = xlsShtNM;
xlswrite(fnameXLS,shtNM,1);

% data
for i = 1:4
    xlswrite(fnameXLS,CELLS{i},i+1);
end



























% display 
%breakpoint at line 1514 peakfit([X;N],centerAt,windowSize,npeaks) --> save('peakfitPlotData','xoffset','baseline','xxx','AA','height')
close all;
runTheRest = 1;
%% extract data
if ~runTheRest
    npeaks = 7;
    % histogram fit window
    centerAt = 2e4;windowSize=4e4;
    peakfit([X;N],centerAt,windowSize,npeaks) % put break point before run and quit and run the rest
    % run --> save('peakfitPlotData','xoffset','baseline','xxx','AA','height');    return
    
end


% clear breakpoint
%%
load rintHist2 % histogram
load peakfitPlotData % fit (xxx,AA,height ...)
recruitmentTraceHistFitFN = 'recruitmentTraceHistFit.tif';



%% fit overlay
figure(11)
hist(rintHist2,100);
hold on;


for m = 1:npeaks
    plot(xxx+xoffset,height(m)*AA(m,:)+baseline,'r','Linewidth',2)  % Plot the individual component peaks in green lines
end

maximize;
    set(gcf,'color','w');  
    set(gca, 'box','off')
grid minor

 
        
%% fit data
figure(12);
[N,X]=hist(rintHist2,100);
[FitResults,FitError]=peakfit([X;N],centerAt,windowSize,npeaks)

figure(11)
positions =FitResults(:,2)
height
set(gca,'Units','normalized')
pos = get(gca, 'Position') % [0.1300, 0.1100, 0.7750, 0.8150] (default)
 
ylim([0 max(N)+10])
for m = 1:npeaks
    x = positions(m);
    x = [x x];
    y = [max(N)+10 max(N)+2]
    xn = (x - min(xlim))/diff(xlim) * pos(3) + pos(1);
    yn = (y - min(ylim))/diff(ylim) * pos(4) + pos(2);
    annotation('textarrow',xn,yn);
end
set(gca, 'FontSize', 24)
       imgFig = getframe(gcf); 
        figCap = imgFig.cdata;
        imwrite(figCap,recruitmentTraceHistFitFN) 

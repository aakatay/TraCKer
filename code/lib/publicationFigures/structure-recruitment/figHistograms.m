% run compileEdgeResults.m after running findEdgeDistanceDistribution.m

isSaveFigs = 0;

cellType = 'HeLa';

fig1FN = sprintf('%s_fig1_edgeDistRatio-Uniform_vs_AreaHist.fig',cellType);
fig2FN = sprintf('%s_fig2_edgeDistRatio-Intensity_vs_AreaHist.fig',cellType);

if isSaveFigs
    figure(101)
    savefig(fig1FN)
    figure(201)
    savefig(fig2FN)
end
    
close all;
openfig(fig1FN); 
openfig(fig2FN); 

imgZFout = [fig1FN(1:end-4) '.tif'];
imgZFout2 = [fig2FN(1:end-4) '.tif'];
delete(imgZFout);
delete(imgZFout2);


%% process border extension
figure(1);% recs
GCA=gca;
GCALine=GCA.Children;
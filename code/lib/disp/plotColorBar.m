function [ticks64, tickPos, ticks, hFig ]  = plotColorBar(varargin)
    
    Cdata = cell2mat(varargin(1));
    nTicks = cell2mat(varargin(2));
    if nargin < 3
        isDisp = 1;
    else
        isDisp = cell2mat(varargin(3));
    end
    MN = min(Cdata(:));
    MX = max(Cdata(:));
    R = MX-MN;
    test = MN:R/999:MX;
    sData = scaleData(test); % scaled data
    sData = uint16(sData);
    MN = 1;
    MX = 64; % size of the colormap
    R = MX-MN;
    tickPos = round(MN:R/(nTicks-1):MX);
    for i = 1:nTicks
        ix = find(sData == tickPos(i));
        IX(i) = round(mean(ix));
    end
    ticks = test(IX);
    
    if ~isDisp, ticks64=0;return; end;
    
    tickPos64 = 1:64;
    for i = 1:64
        ix = find(sData == tickPos64(i));
        IX(i) = round(mean(ix));
    end
    ticks64 = test(IX);
    hFig=figure;
    left=100; bottom=100 ; width=20 ; height=500;
    pos=[left bottom width height];
    axis off
    hC=colorbar([0.1 0.1  0.2  0.8]);
    set(hFig,'OuterPosition',pos) 
    %set(hC,'Ytick',tickPos,'YTicklabel',ticks);
    ticks= round(ticks);
    set(hC,'Ytick',tickPos,'YTicklabel',ticks);
    
    
    
    function sData = scaleData(tData)
        % insert here the scaling equaation of the data


        Cmin = min(tData(tData~=0));
        Cmax = max(tData(:));
        Crange = Cmax - Cmin;
        sData = round(63*(tData-Cmin) / Crange) + 1;
        %sData = tData;
    end
end
    
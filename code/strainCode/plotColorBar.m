function ticks64 = plotColorBar(varargin)

    Cdata = cell2mat(varargin(1));
    nTicks = cell2mat(varargin(2));
    YLabel = [];
    if nargin == 3
        YLabel = cell2mat(varargin(3));
    elseif nargin == 4
        YLabel = cell2mat(varargin(3));
        decimal = cell2mat(varargin(4));
        if decimal == 0, decimal = 1; end
    end
    
    MN = min(Cdata(Cdata~=0));
    MX = max(Cdata(:));
    R = MX-MN;
    tickStep_ = (R/(nTicks));
    dataDigit = floor(log10(tickStep_));
    dataMNdigit = floor(log10(tickStep_));
    tickStep = floor(tickStep_/ power(10,dataDigit)/decimal)*power(10,dataDigit)*decimal;
    MNround = floor(MN/ power(10,dataMNdigit))*power(10,dataMNdigit);
    nTicks = floor(MX/tickStep);
    tickPos_Data = [1:nTicks]*tickStep + MNround;
    colorStep = R/63; % color step
    colorCode = [0:63]*colorStep + MN;
    for i = 1:nTicks
        [~,  tickPos_Color_] = min(abs(colorCode - tickPos_Data(i)));
        tickPos_Color(i) = tickPos_Color_;
    end
    
    fig1=figure;
    left=100; bottom=100 ; width=20 ; height=500;
    pos=[left bottom width height];
    axis off
    hC=colorbar([0.1 0.1  0.2  0.8]);
    set(fig1,'OuterPosition',pos) 
    set(hC,'Ytick',tickPos_Color,'YTicklabel',tickPos_Data);
    set(get(hC,'YLabel'),'String',YLabel)
    
     
    function sData = scaleData(tData)
        % insert here the scaling equaation of the data


        Cmin = min(tData(tData~=0));
        Cmax = max(tData(:));
        Crange = Cmax - Cmin;
        sData = round(63*(tData-Cmin) / Crange) + 1;
        %sData = tData;
    end
end
    
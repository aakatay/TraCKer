function rtSNRvoronoiIMG
% displays SM centers in voronoi cells

    cfg = '..\..\..\cfgRT';
    c_=load(cfg);
    cfg = c_.cfg;
    
    ndigit = cfg.ndigit; % # of digits for sequence number
    label = cfg.label;
    w = cfg.w;
    h = cfg.h;
    mag = cfg.dispMag; % display size
    szXY = [w h];

    [mag, pos, szx, szy ] = calcMaxMag(zeros(szXY),mag);
    szXYmag = [szx, szy];
    szYXmag = fliplr(szXYmag);
    
    
    %% SNR image figure
    tit = 'SNR voronoi image';
    pos(1) = pos(1) - 800;
    pos(1) = pos(1)- 700;        
    figSNRvoronIMG = figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos/2 szXYmag(1) szXYmag(2)]);
    axeSNRvoron = axes('Parent',figSNRvoronIMG,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren','XLim',0.5+[0 szXYmag(1)],'YLim',0.5+[0 szXYmag(2)]);


    % input
    SNRdataFN = 'SNRdata.mat'; % 'ixFrm','XYS'
    
    nlast = 0;
    while (1)
        while (1) % wait for update
            SNRdata_ = dir(SNRdataFN);
            if isempty(SNRdata_), continue; end
            SNRdata=load(SNRdata_.name);
            ixFrm = SNRdata.ixFrm;
            nInput = numel(ixFrm);
            if nlast == nInput
                continue
            else
                nlast = nInput;
                XYS = SNRdata.XYS;
                break; % continue
            end 
        end
        n = nlast;

        CM = gray(256);
        xys = XYS(ixFrm(n-1)+1:ixFrm(n),:);
        if size(xys,1)<5 % write a blank image
            SNRmovVoronoi = zeros(szXY*mag);
            %imwrite(SNRmovVoronoi,SNRmovieVoronoiFN,'WriteMode','append','Compression', 'none') 
        else
            snrDisp = xys(:,3);

            % SNR intensity range scales
            if max(snrDisp)<1.5 
                sscale = 1.5;
                EdgeColorSel = 1; % red
            elseif max(snrDisp)<5
                sscale = 5;
                EdgeColorSel = 2; % green 
            else % max(snrDisp)<10
                sscale = 10;
                EdgeColorSel = 3; % blue
            end
            [SNRmovVoronoi,tPatch(n)] = getVoronoinImg(figSNRvoronIMG,xys(:,1:2),snrDisp/sscale,szXY,mag,CM,EdgeColorSel);
            %imwrite(SNRmovVoronoi,SNRmovieVoronoiFN,'WriteMode','append','Compression', 'none') 
        end
    
    end
    
end
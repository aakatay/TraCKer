function [mnmx] = dispFluorescenceSpectra(varargin)
% called by callFluorescenceSpectra.m
% m : emission
% x : excitation
% FX : filter excitation
% FD : filter dichroic
% FDq : filter dichroic quadview
% FM : filter emission
% xFMq: filter emission quadview

    fc = varargin{1};
    fq = varargin{2};
    fl = varargin{3};
    ls = varargin{4};
    disp = varargin{5};
    isBRT = varargin{6};
    isCy3 = varargin{7};
    isInterp = 0;
    if nargin > 7 % calcFluorescenceSpectra
        mnmx = varargin{8};
        mn = mnmx(1);
        mx = mnmx(2);
        isInterp = 1;
    end
    
    FaceAlphaVal = disp.FaceAlphaVal;
    EdgeAlphaVal = disp.EdgeAlphaVal;
    
    % output
    imgFout = 'fluoSpectra.tif';
    FIGoutFN  = 'fluoSpectra.fig';
    if isBRT
        imgFout = ['BRT-' imgFout];
        FIGoutFN = ['BRT-' FIGoutFN];
    else
        imgFout = ['PRB-' imgFout];
        FIGoutFN = ['PRB-' FIGoutFN];
    end

    if isCy3
        imgFout = ['CY3-' imgFout];
        FIGoutFN = ['CY3-' FIGoutFN];
    else
        imgFout = ['GFP-' imgFout];
        FIGoutFN = ['GFP-' FIGoutFN];
    end 

    %% read input structs
    % (fc) filter cube (Chroma 89401 - ET - DAPI/FITC/TRITC/CY5 Quad)
    fcX = fc.X; % spectra (x-axis)
    fcD = fc.D; % values (y-axis)
    colFC = [0 0 1; 1 0 0; 0 0 0]; % colors
    colFQM = zeros(3,3);
    colFQD = zeros(3,3);
    colFQD(:,1) = 1;
    % (fq)
    fqm = fq.fqm;
    fqd = fq.fqd;
    
    % (fl) fluorophores
    flX = fl.X; % spectra (x-axis)
    flDx = fl.Dx; % excitation values (y-axis)
    flDm = fl.Dm; % emission values (y-axis)
    % (ls) laser sources
    laserX = ls.X;

    %% interpolation
    if isInterp
        xq = mn:0.1:mx; % interp data coors
        xq = xq';
        delete(701)
        figure(702);clf;
    else
        figure(701);clf;
    end
    
    %% display
    maximize
    subplot(3,1,1)
    title('Lasers,Fluos,FilterCube')
    
    % (fc) filter cube
    hold on
    minfcX = 0;
    maxfcX = 10000;
    for i = 1:numel(fcX) 
        if isInterp
            ix =  find((mx >= fcX{i}) .* (fcX{i}>=mn));
            fcD{i} = fcD{i}(ix);
            fcX{i} = fcX{i}(ix);
            fcD{i} = interp1q(fcX{i},fcD{i},xq); %yq = interp1q(x,y,xq);
            fcX{i} = xq;
        else
            minfcX = max([minfcX; min(fcX{i})]);
            maxfcX = min([maxfcX; max(fcX{i})]);
        end
        if mean(fcD{i}) == 1, continue; end % empty
        area(fcX{i},fcD{i},'FaceColor',colFC(i,:),'FaceAlpha',FaceAlphaVal,'EdgeAlpha',EdgeAlphaVal);
    end
    hold off
    
    % (fqd) quadview/CSU
    subplot(3,1,2)
    title('Lasers,Fluos,ExternalDichroic')
    
    minfqdX = 0;
    maxfqdX = 10000;
    if isfield(fqd,'X') % Quadview only
        hold on
        for i = 1:numel(fqd.X) 
            if isInterp
                ix =  find((mx >= fqd.X{i}) .* (fqd.X{i}>=mn));
                fqd.D{i} = fqd.D{i}(ix);
                fqd.X{i} = fqd.X{i}(ix);
                fqd.D{i} = interp1q(fqd.X{i},fqd.D{i},xq); %yq = interp1q(x,y,xq);
                fqd.X{i} = xq;
            else
                minfqdX = max([minfqdX; min(fqd.X{i})]);
                maxfqdX = min([maxfqdX; max(fqd.X{i})]);
            end
            meanF = mean(fqd.X{i}.*fqd.D{i});
            colFQD = spectra([meanF meanF 1]);
            area(fqd.X{i},fqd.D{i},'FaceColor',colFQD,'FaceAlpha',FaceAlphaVal,'EdgeAlpha',EdgeAlphaVal);
        end
        hold off
    end
    
    % (fqm) quadview/CSU
    subplot(3,1,3)
    title('Lasers,Fluos,ExternalEmissionFilter')
    hold on
    minfqmX = 0;
    maxfqmX = 10000;
    for i = 1:numel(fqm.X) 
        if isInterp
            ix =  find((mx >= fqm.X{i}) .* (fqm.X{i}>=mn));
            fqm.D{i} = fqm.D{i}(ix);
            fqm.X{i} = fqm.X{i}(ix);
            fqm.D{i} = interp1q(fqm.X{i},fqm.D{i},xq); %yq = interp1q(x,y,xq);
            fqm.X{i} = xq;
        else
            minfqmX = max([minfqmX; min(fqm.X{i})]);
            maxfqmX = min([maxfqmX; max(fqm.X{i})]);
        end
        meanF = mean(fqm.X{i}.*fqm.D{i});
        colFQM = spectra([meanF meanF 1]);
        area(fqm.X{i},fqm.D{i},'FaceColor',colFQM,'FaceAlpha',FaceAlphaVal,'EdgeAlpha',EdgeAlphaVal);
    end
    hold off
    
    
    for pp = 1:3 % both in two subplots
        subplot(3,1,pp)
        hold on
        % (fl) fluo
        minflX = 0;
        maxflX = 10000;
        fluoN = numel(flX);
        for i = 1:fluoN  
            if isInterp
                ix =  find((mx >= flX{i}) .* (flX{i}>=mn));
                flX{i} = flX{i}(ix);

                flDx{i} = flDx{i}(ix);
                flDx{i} = interp1q(flX{i},flDx{i},xq); %yq = interp1q(x,y,xq);
                flDm{i} = flDm{i}(ix);
                flDm{i} = interp1q(flX{i},flDm{i},xq); %yq = interp1q(x,y,xq);

                flX{i} = xq;
            else
                minflX = max([minflX; min(flX{i})]);
                maxflX = min([maxflX; max(flX{i})]);
            end

            [~,mxIx] = max(flDm{i});
            mxM = flX{i}(mxIx);
            colFl = spectra([mxM mxM 1]);
            plot(flX{i},flDx{i},'Color',[0 0 0],'LineStyle','--','LineWidth',1);
            plot(flX{i},flDx{i},'Color',colFl,'LineStyle','--');
            area(flX{i},flDm{i},'FaceColor',colFl,'FaceAlpha',FaceAlphaVal,'EdgeAlpha',EdgeAlphaVal);

        end

        % (ls) laser
        laserN = numel(laserX);
        for i = 1:laserN
            lx = laserX(i);
            colLaser = spectra([lx lx 1]);
            line([lx lx],[0 1],'Color',colLaser,'LineWidth',2);
        end
        hold off;
    end
    
    
    % calcFluorescenceSpectra
    if isInterp % 2nd call
        mnmx = [];
        calcFluorescenceSpectra
        
        figure(702)
        savefig(FIGoutFN)
        imgFig = getframe(gcf);
        imgOut = imgFig.cdata;
        imwrite(imgOut,imgFout);
    else
        minx = max([minfcX minflX minfqdX]);
        maxx = min([maxfcX maxflX maxfqdX]);
        mnmx = [minx maxx];
    end
    
    
end


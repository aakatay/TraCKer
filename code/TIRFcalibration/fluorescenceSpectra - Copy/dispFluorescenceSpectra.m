function [mnmx] = dispFluorescenceSpectra(varargin)

    fc = varargin{1};
    fl = varargin{2};
    ls = varargin{3};
    disp = varargin{4};
    isInterp = 0;
    if nargin > 4
        mnmx = varargin{5};
        mn = mnmx(1);
        mx = mnmx(2);
        isInterp = 1;
    end
    
    FaceAlphaVal = disp.FaceAlphaVal;
    EdgeAlphaVal = disp.EdgeAlphaVal;

    %% INPUTs
    % filter cube (Chroma 89401 - ET - DAPI/FITC/TRITC/CY5 Quad)
    fcX = fc.X; % spectra (x-axis)
    fcD = fc.D; % values (y-axis)
    colFC = [0 0 1; 1 0 0; 0 0 0]; % colors
    % fluorophores
    flX = fl.X; % spectra (x-axis)
    flDx = fl.Dx; % excitation values (y-axis)
    flDm = fl.Dm; % emission values (y-axis)
    % laser sources
    laserX = ls.X;


    %% interpolation
    if isInterp
        xq = mn:0.1:mx; % interp data coors
        xq = xq';
        figure(701);clf;
    else
        figure(702);clf;
    end
    
    %% display
    maximize

    hold on
    minfcX = 0;
    maxfcX = 10000;
    for i = 1:numel(fcX)
        if isInterp
            ix =  find((mx >= fcX{i}) .* (fcX{i}>=mn));
            fcDy = fcD{i}(ix);
            x = fcX{i};
            y = fcD{i};
            yq = interp1q(x,y,xq);
            fcX{i} = xq;
            fcD{i} = yq;
            
        else
            minfcX = max([minfcX; min(fcX{i})]);
            maxfcX = min([maxfcX; max(fcX{i})]);
        end
        area(fcX{i},fcD{i},'FaceColor',colFC(i,:),'FaceAlpha',FaceAlphaVal,'EdgeAlpha',EdgeAlphaVal);
    end
    
    minflX = 0;
    maxflX = 10000;
    for i = 1:numel(flX)
        [~,mxIx] = max(flDm{i});
        mx = flX{i}(mxIx);
        colFl = spectra([mx mx 1]);
        plot(flX{i},flDx{i},'Color',colFl,'LineStyle','--');
        area(flX{i},flDm{i},'FaceColor',colFl,'FaceAlpha',FaceAlphaVal,'EdgeAlpha',EdgeAlphaVal);
        
        minflX = max([minflX; min(flX{i})]);
        maxflX = min([maxflX; max(flX{i})]);
    end
    
    for i = 1:numel(laserX)
        lx = laserX(i);
        colLaser = spectra([lx lx 1]);
        line([lx lx],[0 1],'Color',colLaser,'LineWidth',2);
    end
    hold off;
    
    if isInterp
        mnmx = [];
    else
        minx = max([minfcX minflX]);
        maxx = min([maxfcX maxflX]);
        mnmx = [minx maxx];
    end
    cxcc=2;
    
end


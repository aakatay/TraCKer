function [IMG,tPatch] = getVoronoinImg(figSNRvoronIMG,x,cellColor,szXY,mag,CM,EdgeColorSel)
% displays scattering centers in voronoi cells
% test : getVoronoinImg(gallery('uniformdata',[10 2],5)*50,[50 50],8)
tPatch = [];
    isSymmetricCorrection = 1;
    isDisp = 0; 
    EdgeColor = [1 0 0; 0 1 0; 0 0 1]; % red green blue
    
    % check and correct direction of input
    if size(x,1)<size(x,2) % rows
        x = x';
    end
    if size(szXY,1)>size(szXY,2) % columns
        szXY = szXY';
    end
    
    % normalize the coors
    x = x ./ repmat(szXY,size(x,1),1);

    %% first call
    if isDisp, figure(1); end
    [v,c] = voronoin(x); 
    ix = [];
    for i = 1:length(c) 
        if all(c{i}~=1)   % If at least one of the indices is 1, 
                          % then it is an open region and we can't 
                          % patch that.
            if isDisp
                p=patch(v(c{i},1),v(c{i},2),i); % use color i.
            end
        else
            ix = [ix i]; % indices of non-bounded cells
        end
    end
    if isDisp
        hold on;
        scatter(x(:,1),x(:,2),'.')
        scatter(v(:,1),v(:,2),'.','r')
        hold off
    end

    if isSymmetricCorrection
    %% find empty quarter
        qix = (x(ix,1)>0.5)*10 +  (x(ix,2)>0.5)*100; % quarter index
        % 0: bottom left
        % 10: bottom right
        % 100: top left
        % 110: top right
        qixs = [0 10 100 110];
        corners = [0 0; 1 0; 0 1; 1 1];
        qem = [];
        for i = 1:4
            if ~ismember(qixs(i),qix)
                qem = [qem i]; % missing quarters
            end
        end


    %% find one center for each corner    
        xr = x; % remaining centers
        ixr = ix; % remaining non-bounded cells

        % 1: assign centers to corners where quarters are empty
        u = x(ixr,1); % x coordinate
        v = x(ixr,2); % y coordinate
        xcrnr = nan(1,4);
        ixCorner = [];
        for i = 1:numel(qem)
            cx = corners(qem(i),1);
            cy = corners(qem(i),2);
            % find closest to corner
            d = sqrt((u-cx).^2 + (v-cy).^2); % distance
            [~,ixmin] = min(d);
            xr(ix(ixmin),:) = nan;
            xcrnr(qem(i)) = ixr(ixmin);
            ixCorner = [ixCorner ixr(ixmin)];
            ixr(ixmin) = [];
            u(ixmin) = [];
            v(ixmin) = [];
        end


        % 2: assign centers to remaining corners
        u = xr(ixr,1); % x coordinate
        v = xr(ixr,2); % y coordinate
        xcrnrRem = find(isnan(xcrnr)); % remaining corners
        for i = 1:numel(xcrnrRem)
            cx = corners(xcrnrRem(i),1);
            cy = corners(xcrnrRem(i),2);
            % find closest to corner
            d = sqrt((u-cx).^2 + (v-cy).^2); % distance
            [~,ixmin] = min(d);
            xr(ixr(ixmin),:) = nan;
            xcrnr(xcrnrRem(i)) = ixr(ixmin);
            ixCorner = [ixCorner ixr(ixmin)];
            ixr(ixmin) = [];
            u(ixmin) = [];
            v(ixmin) = [];
        end

    %% generate symmetrical centers outside the image
        % 1: 4 new centers symmetrical to corners
        x2 = x;
        for i = 1:4 % each corner
            xnew = 2*corners(i,:) - x(xcrnr(i),:);
            x2 = [x2; xnew];
        end

        % 2: the rest symmetrical to edges
        x3 = x2;
        % edges
        %1: x=0;
        %2: x=1;
        %3: y=0;
        %4: y=1;
        sx = [-1 -1 1 1];
        sy = [1 1 -1 -1];
        dx = [0 2 0 0];
        dy = [0 0 0 2];
        u = xr(ixr,1); % x-coor
        v = xr(ixr,2); % y-coor

        [~,ixEdge] = min( [u (1-u) v (1-v)] );
        unew = dx(ixEdge)+ u.*sx(ixEdge);
        vnew = dy(ixEdge)+ v.*sy(ixEdge);
        xnew = [unew' vnew'];
        if isempty(xnew), xnew=[];end
        x3 = [x3; xnew];
    else
        x3 = x;
    end



%    t1 = toc;
    %% second call
    m = szXY(1)*mag;
    n = szXY(2)*mag;
    tit = 'getVoronoinImg';
    pos2 = [300 300];
    %figImg = figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',CM,'Position',[pos2/2 m n]);
    %axeImg = axes('Parent',figImg,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','XLim',0.5+[0 m],'YLim',0.5+[0 n]);
    
    xlim3 = [0 1];
    ylim3 = [0 1];
    xlim(xlim3)
    ylim(ylim3)
    [v,c] = voronoin(x3); 
    nc = size(c,1);
    cellColor2 = zeros(nc,1);
    cellColor2(1:numel(cellColor)) = cellColor;
    for i = 1:length(c) 
        if all(c{i}~=1)   % If at least one of the indices is 1, 
                          % then it is an open region and we can't 
                          % patch that.
            p = patch(v(c{i},1),v(c{i},2),[1 1 1]*cellColor2(i),'EdgeColor',EdgeColor(EdgeColorSel,:)); % use color i.
            %p = patch(v(c{i},1),v(c{i},2),[1 1 1]*cellColor2(i)); % use color i.
            cccc = 33;
        end
    end
    hold on;
    scatter(x3(:,1),x3(:,2),'.','r')
    %scatter(v(:,1),v(:,2),'.','r')
    hold off
    
    %% 
    figure(figSNRvoronIMG);
    imgFig = getframe(gcf); 
    IMG = imgFig.cdata;
    IMG = flipud(IMG);
    %close(figImg);
%    t2 = toc;
%    tPatch = t2-t1;


    %% show extra generated centers
    if isDisp 
        xlim2 = xlim;
        ylim2 = ylim;
        figure(1)
        xlim(xlim2)
        ylim(ylim2)
        
        figure(3);
        hold on;
        n1 = size(x,1);
        n2 = size(x2,1);
        scatter(x(:,1),x(:,2),'.','b')
        scatter(x3(n1+1:n2,1),x3(n1+1:n2,2),'.','g')
        scatter(x3(n2+1:end,1),x3(n2+1:end,2),'.','r')
        scatter(x3(ixCorner,1),x3(ixCorner,2),'.','m')

        xlim(xlim2)
        ylim(ylim2)
        h = imrect(gca, [0 0 1 1]);
        hold off
    end



    %% print image
    
    % image1: overlay     
end






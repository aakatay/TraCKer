    function dispBoundaries(R,Bt,Lt,mag,ps,isSelSt)
    if isempty(ps), ps = nan(1,3); end
    imagesc(R); axis image;
    CM = colormap;
    CM(1,:)=0.1;
    colormap(CM)
    
    hold on     
    for k = 1:length(Bt) % each tight section
        boundaryT = Bt{k};
        
        isdrAmin = 0;
        genStructRec; %R,Lt -> dr,Nr generates structures  to filter out skinnies
        
        % color processed structs
        psix = find(ps(:,1) == k);
        cB = 0;
        if ~isempty(psix) 
            if isSelSt, cL = [0.5 0.5 0.5]; else, cL = 'g'; end; % default color is gray or green
            cC = ps(psix,2); % color code
            if cC>=20, cB=1; cC=cC-20; end; % bold
            switch cC % color code -> color letter
                case 1
                    cL = 'k'; 
                case 2
                    cL = 'w';
                case 3
                    cL = 'b';
                case 4
                    cL = 'c';
                case 5
                    cL = 'm'; 
                case 6
                    cL = 'y';
                case 11
                    cL = 'g';
            end
            %black white blue cyan yellow magenta
        else
            cL = 'g'; % not processed
        end
        bw=1;
        [Nr,~]=histcounts(drmin,'BinWidth',bw);
        skCoeff = 0.5; % rate of edge points
        if size(boundaryT,1)==2 % point
            scatter(boundaryT(:,2), boundaryT(:,1),'*','r')
        elseif 0 && Nr(1)/sum(Nr)>=skCoeff % skinny structure
            plot(boundaryT(:,2), boundaryT(:,1), 'r', 'LineWidth', 1+cB)
        else
            plot(boundaryT(:,2), boundaryT(:,1),'color', cL, 'LineWidth', 1+cB)
        end
    end
    hold off
    set(gcf,'units','pixels','Position',[120,120,size(R,2)*mag,size(R,1)*mag]); 
    set(gca,'units','pixels','Position',[0,0,size(R,2)*mag,size(R,1)*mag]);
    
    stix = 1;
    if isSelSt & size(ps,1)>10 % display structure numbers
        figure(111)
        imagesc(R.*0); axis image;
        for i = 1:size(ps,1)
            if ps(i,2)==1 | ps(i,2)==11
                [yy,xx] = find(Lt(:,:,i)>0);
                x=mean(xx)-10; y=mean(yy)-10;
                if x <= 0, x=1; end
                if y <= 0, y=1; end
                text(x,y,sprintf('%i',stix));
                stix = stix+1;
            else
                continue;
            end
        end
        set(gcf,'units','pixels','Position',[120,120,size(R,2)*mag,size(R,1)*mag]); 
        set(gca,'units','pixels','Position',[0,0,size(R,2)*mag,size(R,1)*mag]);
    end
    if isempty(findobj('type','figure','Number',888))
        paintStructs_w_EdgeDistance_dontPrint=1;
        figure(888)
        paintStructs_w_EdgeDistance;
        set(gcf,'units','pixels','Position',[120,120,size(R,2)*mag+100,size(R,1)*mag]); 
        set(gca,'units','pixels','Position',[0,0,size(R,2)*mag,size(R,1)*mag]);
    end
end
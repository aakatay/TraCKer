function updCoeff(hObject,event,hFig,hTextCoeff,hscat,hscat2,hscat3,CoeffFit,sizeImg,nRow,IMGtiled,IMGfiltTiled,filtMax,isDebugFilt,frstFrm,nTile,nFrTile)


    slideVal = get(hObject,'Value');
    Coeff = slideVal*filtMax;
    framesVec = frstFrm + [(nTile-1)*nFrTile : nTile*nFrTile-1];
    for j = framesVec
        k = j - frstFrm -(nTile-1)*nFrTile + 1; % frame number in the tile
        CoeffFitTile_(:,:,k) = ones(sizeImg)*CoeffFit(j);
    end
    % find [x y] after coeff update
    CoeffFitTiled = tileFrames(CoeffFitTile_,nRow);
    CoeffTiled = CoeffFitTiled*Coeff;
    
    
    % DETECTION THRESHOLD
    IMGfilt=IMGfiltTiled; IMG=IMGtiled; CoeffThresh=CoeffTiled;
    detectThreshold; % (IMGfilt,IMG,CoeffThresh) --> (BW,din)
        
    [y,x,v]=find(BW==1); 
    if numel(BW(:)) == sum(BW(:))
        x=[]; y=[];
    end;
        
    set(hscat,'Xdata',x,'Ydata',y);
    
    if isDebugFilt
        set(hscat2,'Xdata',x,'Ydata',y);
        set(hscat3,'Xdata',x,'Ydata',y);
    end
    figure(hFig);
    set(hTextCoeff,'String',sprintf('coefficient: %i',round(Coeff)));
    
end
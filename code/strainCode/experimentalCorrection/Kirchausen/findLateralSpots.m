function BasoLateral = findLateralSpots(mxPj,MxZeInd,nSpot)
% returns 
% mxPj,MxZeInd outputs of the max proj. of the ventral section of the cell

    [mxMxPj mxMxPjInd] = max(mxPj(:));
        
    %define spot window for filtering out spots
    sZ = 2; %spot size (in px)
    sP = ones(sZ*2+1,sZ*2+1); %spot window
    
    [sx sy] = size(mxPj);
    
    mxPjFlt = 0.*mxPj(:);
    mxPjSpt = 0.*mxPj(:);
    mxPj2   = mxPj;
    for iN = 1:nSpot % # of spots analyzed
        [mxMxPj mxMxPjInd] = max(mxPj2(:));
        [srtMxMxPj srtMxMxPjInd]= sort(mxMxPj);
        %% save X Y Z coord of the spot
        [Ix(iN) Iy(iN)] = linInd2xyInd(mxMxPjInd(srtMxMxPjInd(1)),sy,sx); %converts linear indice to 2D xy indice 
        Iz(iN) = MxZeInd(Ix(iN),Iy(iN));
        
        %% filter out the max spot
        mxPjFlt(mxMxPjInd(srtMxMxPjInd(1))) = 1; % mark the center of the spot
        mxPjFlt2Dind    = reshape(mxPjFlt,size(mxPj,1),size(mxPj,2));        
        mxPjFlt2D       = conv2(mxPjFlt2Dind,sP);
        mxPjFlt2D       = mxPjFlt2D(sZ+1:end-sZ,sZ+1:end-sZ); % resize by cropping after convolution
        mxPjFlt2D       = 1-mxPjFlt2D; % invert for substraction
        %mxPjSpt         = mxPjSpt + mxPjFlt2D; % maximum proj spots
        mxPj2           = mxPjFlt2D.*mxPj2; % remove the spot       
        %figure(3); imagesc(mxPj2);
        %figure(4); scatter3(Iy',Ix',Iz');
        %display('a')
        
    end
    figure(3); imagesc(mxPj2)
    BasoLateral = [Ix' Iy' Iz'];
    figure(4); scatter3(Iy',Ix',Iz');
end
    
    
function CM = genColorMap(colorMap,Frames)
% generates colormap for coding the time info to color 

    %upsampleRate = 10;
    N = Frames+1; % # of time frames
    nRep = 1; % number of repitiions
    szCM = N /nRep;
    CM0 = colormap(colorMap);
    szCM0 = size(CM0,1);
    s = (size(CM0,1)-1)/szCM;
    
    vv = 1-s:s:szCM0;
    vv = vv(1:szCM);
    CM = interp1(1:szCM0,CM0(:,1),vv)';
    CM(:,2) = interp1(1:szCM0,CM0(:,2),vv);
    CM(:,3) = interp1(1:szCM0,CM0(:,3),vv);
    CM = CM(2:end,:);
    return;
    
    
    CM2 = flipud(CM);
    CM3 = [CM' CM2']';
    CM4 = repmat(CM3,nRep,1);
    close(gcf)
end
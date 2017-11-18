
    figure(11)
    subplot(2,4,1);
    hist(trInf(:,12),100); title('max intensity')
    subplot(2,4,2);
    hist(trInf(:,6),100); title('mean intensity')
    subplot(2,4,3);
    hist(trInf(:,11),100); title('sum intensity')
    subplot(2,4,4);
    scatter(trInf(:,12),trInf(:,7),'.'); title('mean spread vs max intensity')
    subplot(2,4,5);
    scatter(trInf(:,12),trInf(:,11),'.'); title('sum intensity vs max intensity')
    subplot(2,4,6);
    scatter(trInf(:,12),trInf(:,6),'.'); title('mean intensity vs max intensity')
    subplot(2,4,7);
    scatter(trInf(:,8),trInf(:,7),'.'); title('mean spread vs max spread')
    subplot(2,4,8);
    scatter(trInf(:,8),trInf(:,2),'.'); title('trace length vs max spread')
    % write figure
    maximize
    imgFig = getframe(gcf);
    tpImg = imgFig.cdata; 
    tpFN = sprintf('trace populations.tif');
    imwrite(uint16(tpImg),tpFN)

    %% set2
    [N,X] = hist(trInf(:,1),100);

    figure(10)
    subplot(2,4,1);
    plot(X,smooth(N)); title('number of traces vs trace time')
    subplot(2,4,2);
    scatter(trInf(:,7),trInf(:,11),'.'); title('sum intensity vs mean spread')
    subplot(2,4,3);
    scatter(trInf(:,8),trInf(:,11),'.'); title('sum intensity vs max spread')
    subplot(2,4,4);
    scatter(trInf(:,1),trInf(:,12),'.'); title('max intensity vs trace time')
    subplot(2,4,5);
    scatter(trInf(:,1),trInf(:,6),'.'); title('mean intensity vs trace time')
    subplot(2,4,6);
    scatter(trInf(:,1),trInf(:,11),'.'); title('sum intensity vs trace time')
    subplot(2,4,7);
    scatter(trInf(:,8),trInf(:,7),'.'); title('mean spread vs max spread')
    subplot(2,4,8);
    scatter(trInf(:,8),trInf(:,2),'.'); title('trace length vs max spread')
    % write figure
    maximize
    imgFig = getframe(gcf);
    tpImg = imgFig.cdata; 
    tpFN = sprintf('trace populations2.tif');
    imwrite(uint16(tpImg),tpFN)
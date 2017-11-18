function plotBoundaries(B,A,isSc)
% display boundaries in overlay image
    if isSc
        imagesc(A);
    else
        image(A);
    end
    hold on
    [szy,szx] = size(A);
    for i = 1:numel(B)
        c=cell2mat(B(i));
        plot(c(:,2),c(:,1),'LineWidth',1,'Color','r');
    end
    hold off
    axis image
    xlim([0.5 szx+0.5])
    ylim([0.5 szy+0.5])
        
       
end
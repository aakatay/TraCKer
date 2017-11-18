% plots histogram with logarithmic count axis
function [histcount,centers] = loghist(img,n,dbg)
% n : number of bins


    [histcount,edges] = histcounts(img(:),n);
    edgesDiff = edges(2) - edges(1);
    centers = edges + edgesDiff/2;
    centers(end) = [];
    if dbg,    semilogy(centers,histcount+1); end;
    %plot(centers,histcount+1)
end
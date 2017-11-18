function mxDist = findMaxDist(TRACEX,TRACEY)
%   finds the max distance betw. xy locations of the trace
    TRACEX(TRACEX==0) = nan;
    TRACEY(TRACEY==0) = nan;
    n = numel(TRACEX);
    tracey1 = repmat(TRACEY,1,n);
    tracey2 = reshape(tracey1,n,n);
    tracey2 = tracey2';
    tracey2 = tracey2(:)';
    tracex1 = repmat(TRACEX,1,n);
    tracex2 = reshape(tracex1,n,n);
    tracex2 = tracex2';
    tracex2 = tracex2(:)';
    
    
    
    mxDist = max(sqrt((tracex1-tracex2).^2+(tracey1-tracey2).^2)); % MAX deviation
    
end
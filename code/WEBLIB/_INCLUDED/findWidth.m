function FWHM = findWidth(fx)
%findFWHM: Finds the full width half max (FWHM) of a function
% nmx : number of max values to find a max value for the function
%	FWHM = findFWHM(x, fx);

    fx = double(fx);
    
    fxDiff = fx(1:end-1)-fx(2:end);
    %figure; plot(fxDiff)
    [~,mnix] = min(fxDiff);
    [~,mxix] = max(fxDiff);
    mnix = mnix-0.5;
    mxix = mxix+0.5;
    
    if 0 

        ind = find(fx>=m/2);	%	Find indicies where I>=max(I)/2
        nl = min(ind);			%	Leftmost index
        nr = max(ind);			%	Rightmost index

        %	Linear interpolate x positions
        xl = (x(nl)-x(nl-1))  *  (m/2-fx(nl-1))/(fx(nl)-fx(nl-1)) + x(nl-1);
        xr = (x(nr-1)-x(nr))  *  (m/2-fx(nr))/(fx(nr-1)-fx(nr))  + x(nr);
    end

    %	Get FWHM
    FWHM = mxix-mnix;
end

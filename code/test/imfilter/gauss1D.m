function h = gauss1D(sz,sgm)
% gaussian 1D
% sz: (size) has to be odd
% sgm: (sigma) std of the gaussian func

     siz   = (sz-1)/2;
     std   = sgm;
     x = -siz:siz;
       
     arg   = -(x.*x)/(2*std*std);

     h     = exp(arg);
     h(h<eps*max(h(:))) = 0;

     sumh = sum(h(:));
     if sumh ~= 0,
       h  = h/sumh;
     end;
end
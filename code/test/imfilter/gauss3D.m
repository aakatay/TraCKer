function h = gauss3D(szXY,sgmXY,szZ,sgmZ,x,y,z)
% gaussian 3D
% sz: (size) actual size = 2*sz+!
% sgm: (sigma) std of the gaussian func
% x

     x = -szXY:szXY;
     z = -szZ:szZ;
     
     [a b c]= ndgrid(x,x,z);
     
     arg   = -(a.*a + b.*b)/(2*sgmXY*sgmXY)-c.*c/(2*sgmZ*sgmZ);

     h     = exp(arg);
     h(h<eps*max(h(:))) = 0;

     sumh = sum(h(:));
     if sumh ~= 0,
       h  = h/sumh;
     end;
end


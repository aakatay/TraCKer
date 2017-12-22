    [X,Y,V] = peaks(10);
    [Xq,Yq] = meshgrid(-6.5:6.5,-3:3);
    Vq = interp2(X,Y,V,Xq,Yq);
    imagesc(Vq)
    figure
    imagesc(V)
    
    
    interp1q([1:5]',[1:5]',[1.1 2.4 3.3 4.4 4.9]')
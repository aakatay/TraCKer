    [X,Y,V] = peaks(10);
    [Xq,Yq] = meshgrid(-6.5:6.5,-3:3);
    Vq = interp2(X,Y,V,Xq,Yq);
    imagesc(Vq)
    figure
    imagesc(V)
function gammaFilter

    N = 100;
    x = 1:N;
    x= abs(x-N/2);
    global X Y;
    [X,Y] = meshgrid(x,x);
    hImg = imagesc(X.^2+Y.^2);
    h = uicontrol('style','slider','units','pixel','position',[20 20 300 20]);
    addlistener(h,'ActionEvent',@(hObject, event) makeplot(hObject, event,x,hImg));




end
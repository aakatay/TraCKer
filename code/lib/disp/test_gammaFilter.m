
function gammaFilter

    fname = 'dataNorm.tif';

    imageInfo=imfinfo(fname);
    numFrames=length(imageInfo);
    %Size the movie matrix
    imSize=[imageInfo(1).Height,imageInfo(1).Width,numFrames];

    frstFrame =1;
    for i = frstFrame : numFrames
        A(:,:,i-frstFrame+1) = double(imread(fname,i));
    end
    numFrames = numFrames - frstFrame+1;

    A2 = A(:,:,1);
    % Adiff = cat(3,A2,ones(imSize(1),imSize(2))) - cat(3,ones(imSize(1),imSize(2)),A2);
    % Adiff = Adiff(:,:,2:numFrames-1);
    n = 1;
    mxA = max(A(:)); thr =0;
    hImg = imagesc(A2.^n-mxA.*thr);
    h = uicontrol('style','slider','units','pixel','position',[20 20 300 20]); % gamm
    h2 = uicontrol('style','slider','units','pixel','position',[330 20 300 20]); % threshold
    addlistener(h,'ActionEvent',@(hObject, event) makeplot(hObject, event,A2,hImg));
    addlistener(h2,'ActionEvent',@(hObject2, event) makeplot2(hObject2, event,A2,hImg));


    function makeplot(hObject,event,x,hImg)
        n = get(hObject,'Value');
        n = (10*n + 1);
        
        

        A3 = A2.^n-mxA.*thr;
        %Adiff = cat(3,A3,ones(imSize(1),imSize(2))) - cat(3,ones(imSize(1),imSize(2)),A3);
        %Adiff = Adiff(:,:,2:numFrames-1);

        set(hImg,'Cdata',A3);
        drawnow;
    end
    function makeplot2(hObject,event,x,hImg)
        thr = get(hObject,'Value');
        n = thr;
        A3 = A2.^n-mxA.*thr;
        %Adiff = cat(3,A3,ones(imSize(1),imSize(2))) - cat(3,ones(imSize(1),imSize(2)),A3);
        %Adiff = Adiff(:,:,2:numFrames-1);

        set(hImg,'Cdata',A3);
        drawnow;
    end
    function makeplot233(hObject,event,x,hImg)
        thr = get(hObject,'Value');
        disp(thr);
        n = thr;
        A3 = A2.^n-mxA.*thr;
        %Adiff = cat(3,A3,ones(imSize(1),imSize(2))) - cat(3,ones(imSize(1),imSize(2)),A3);
        %Adiff = Adiff(:,:,2:numFrames-1);

        set(hImg,'Cdata',A3);
        drawnow;
    end
end
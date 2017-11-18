function Jcv = tileFrames(JF,m)
% tiles the frames in a single image
% m : number of rows

    nFr = size(JF,3); % number of frames
    n =  nFr/m;
    %n=n-1;
    if nFr ~= 1, JF = padarray(JF,[1 1]); end;
    Jcv = [];
    for i = 1:m % rows
        Jch = JF(:,:,(i-1)*n+1);
        for j = (i-1)*n+2:i*n
            Jch = horzcat(Jch,JF(:,:,j));
        end
        Jcv = vertcat(Jcv,Jch);
    end
end
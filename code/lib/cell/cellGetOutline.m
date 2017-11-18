
clear
isViewOnly =0; isEdit = 1;
% marks the outline of the cell through time
fname = dir('maxProjDeep*.*')
CM = 'hot';
%CM = 'autumn';
%CM = 'summer';

%fname = dir('maxProj*.*')
imageInfo=imfinfo(fname.name);
numFrames=length(imageInfo);
%Size the movie matrix
imSize=[imageInfo(1).Height,imageInfo(1).Width,numFrames];

    for z=1:numFrames % # of frames
        I(:,:,z) = imread(fname.name,z);
        % normalization to correct for bleach
        Itot(z) = sum(sum(I(:,:,z)));
        %Itot(z) = 
        if z == 1 
            ItotFirst = Itot;
            Itot = 1;
        else
            Itot(end) = Itot(end)/ItotFirst;
        end
        I(:,:,z) = I(:,:,z)./Itot(end);
    end
    
    % normalize the whole stack
    Imax = double(max(I(:)));
    I = uint8(double(I)*256 / Imax);

    if isViewOnly || isEdit
        load('outlineCoor')
    else
        pos = [];
    end
%%
    q = 0; zIx =1; posIx(1)=1;posIxIx =2;isFirst=1; isPosUpd=0;
    szX = [1 imSize(2)];
    szY = [1 imSize(2)];
    while q == 0
        if zIx < 1, zIx = 1; end;
        if zIx > numFrames, zIx = numFrames; end;
        if zIx == 1 && isFirst
            image(I(:,:,zIx));axis equal;colormap(CM)
            set(gcf,'Color',[0.8 0.8 0.8]);
            if ~isempty(find(posIx==zIx))
                set(gcf,'Color',[1 1 1]);
            end
            zIxPrev = zIx;
            if isempty('pos')
                hPol = impoly();
            else
                hPol = impoly(gca,pos(:,:,zIx));
            end
            isFirst = 0;
        else
            if size(pos,3)<zIx || isPosUpd
                pos(:,:,zIxPrev) = getPosition(hPol);   % write coord.
            end
            if isPosUpd % recalculate positions
                [posIx, pos] = cellGetOutline_calcPos(posIx,pos);
                isPosUpd = 0;
            end
            image(I(:,:,zIx));axis equal;colormap(CM)
            set(gcf,'Color',[0.8 0.8 0.8]);
            if ~isempty(find(posIx==zIx))
                set(gcf,'Color',[1 1 1]);
            end
            if size(pos,3)<zIx
                hPol = impoly(gca,pos(:,:,zIxPrev));
            else
                hPol = impoly(gca,pos(:,:,zIx));
            end
            zIxPrev = zIx;
        end
        title(sprintf('frame#: %i',zIx));
        h = figure(1);
        btn = 0;
        while btn == 0
            btn = waitforbuttonpress;
            k = get(h,'CurrentCharacter');
        end
        switch lower(k)
            case 'a'
                zIx = zIx - 1;
            case 'd'
                zIx = zIx + 1;
            case 's'
                pIx = find(posIx==zIx);
                posIxIx = posIxIx - 1;
                posIx = [ posIx(1:pIx-1) posIx(pIx+1:end) ];
            case 'q'
                zIx = zIx - 5;
            case 'e'
                zIx = zIx + 5;
            case 'w'
                zIx = zIx + 15;
            case 'x'
                wait(hPol);
                posIx(posIxIx) = zIx; posIx = sort(posIx);
                posIxIx = posIxIx + 1;
                isPosUpd = 1; % position updated
            case 's' % select spots
        end
        if zIx > numFrames
            isReview = input('review? (y/n)','s');
            if isReview == 'n'
                q = 1;
            end
        end
        if q == 1
            if zIx > numFrames, zIx = numFrames; end;
            pos(:,:,zIx) = getPosition(hPol);
            BW(:,:,zIx) = roipoly(szX,szY,I(:,:,zIxPrev),pos(:,1),pos(:,2));
        else
            ;
        end
    end
    % 
    [posIx, pos] = cellGetOutline_calcPos(posIx,pos);

    Nfr = size(pos,3); % # frames
    for i = 1:Nfr
        BW(:,:,i) = roipoly(szX,szY,I(:,:,i),pos(:,1,i),pos(:,2,i));
        imagesc(BW(:,:,i)); title(sprintf('frame#: %i',i));
        pause(0.01)
    end
    save('outlineCoor','pos','posIx');
    save('outlineMask','BW');
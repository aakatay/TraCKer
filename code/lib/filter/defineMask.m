        mag = 2;
        A=imread(fname);
        figure(41)
        imagesc(A)
        set(gcf,'units','pixels','Position',[720,-120,size(A,2)*mag,size(A,1)*mag]); 
        set(gca,'units','pixels','Position',[0,0,size(A,2)*mag,size(A,1)*mag]);
        % draw a poly for separating the image two pieces
        hp = impoly(gca);
        pos = getPosition(hp);

        pos1 = zeros(size(pos,1)+2,2); % crop 1
        pos1(2:end-1,:) = pos;
        pos1(1,:) = [0,0];
        pos1(end,:) = [0,size(A,1)];
        mask = roipoly(A,pos1(:,1),pos1(:,2));
        figure(42)
        imagesc(mask)
        save(fnameMask,'pos1','mask');
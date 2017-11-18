%% 1- select a rectangle to remove outside dots (stay in image size)
%% 2- to quit select a rect towards outside on upper left of the image 
%% 3- it will find edge and distance to the edge with reference and plot and write to xls
clear
close all
frec=rdir('rec*');
%isGenStructRec = 1; % generates structures to compare
%%
for j = 1:numel(frec)
    close all;
    fname = frec(j).name;
    fnameP = [ 'rcP' fname(4:end)];
    fnameX = [ 'rcX' fname(4:end)];
    fnameB = [ 'bounds' fname(4:end)]; % boundaries
    fnameD = [ 'edgeDistHist' fname(4:end)]; % distribution
    fnameA = [ 'ST' fname(4:end)]; % distribution
    fnameMAT = [ 'edgeDistData' fname(4:end-3) 'mat']; % distribution
    cnvSz = 5;
    convWin = ones(cnvSz);
    if ~exist(fnameP)  % if outlier recruitments are not removed
        R = double(imread(fname));
        R0 = R;
        while 1 % remove outside data
            Rc = conv2(R,convWin,'same');
            figure(1);imagesc(R); axis image;
            figure(2);imagesc(Rc); axis image;

            [B,L]=bwboundaries(Rc,'noholes');
            figure(3); imagesc(L); axis image;
            hold on 
            [Y,X]=find(R>0);
            scatter(X,Y,'.','r');
            hold off
            h=imrect(gca);
            p = round(getPosition(h));
            if p(1)<=0 || p(2)<=0, break; end  % BREAK
            p4=p(2)+p(4);
            p3=p(1)+p(3);
            if p(2)+p(4)> size(R,1), p4=size(R,1); end
            if p(1)+p(3)> size(R,2), p3=size(R,2); end
            R(p(2):p4,p(1):p3) = 0;
            %get(gcf,'CurrentCharacter')
        end
        imwrite(uint16(R),fnameP);
        Rdiff = R0-R;
        if sum(Rdiff(:)), imwrite(uint16(Rdiff),fnameX); end;
    elseif exist(fnameD) % already done
        continue;
    else
        R = double(imread(fnameP));
            Rc = conv2(R,convWin,'same');

            [B,L]=bwboundaries(Rc,'noholes');
    end
    ix_=strfind(fname,'_');
    %%
    figure(1); clf; % edge selection overlay with rec.
    imagesc(R); axis image;
    hold on 
    ilast = 0;
    ilastr = 0;
    ilastrA = 0;
    d=[];dr=[];drA=[];Bt={};s=1;Lt=[];
    if numel(B)>1, disp('multiple sections'); end;
    hw = waitbar(0,'sectioning');
    for k = 1:length(B) % all sections
        boundary = B{k};
        %plot(boundary(:,2), boundary(:,1), 'k', 'LineWidth', 2)
        %%
        %find each section area
        rb = zeros(size(R));
        for b = 1:size(boundary,1) % draw the boundary on image
            rb(boundary(b,1),boundary(b,2))=1;
        end
        wb =double(watershed(rb,4)); wb=wb-1;wb(wb~=0)=1; 

        % find tighest boundary
        rbThick = conv2(rb,ones(3),'same');
        rbThick(rb>0)=0; rbThick0 = rbThick;
        rbThick(wb==0)=0;    
        wbThick =double(watershed(rbThick,8));
        tBoundsImg = zeros(size(wbThick)); % tight boundary
        tBoundsImg(wbThick>=2)=1; tBoundsImg0 = tBoundsImg;
        tBoundsImg(rbThick>0)=0;
        [Bt0,Lt0]=bwboundaries(tBoundsImg,'noholes'); % tight boundaries
        %figure(11); imagesc(Lt0); axis image
        %figure(12); imagesc(Lt); axis image
        if numel(Bt0)>0 % multiple sub-sections
            if sum(Lt0(:)) < 50, continue; end
            Bt = [Bt;Bt0];
            for i = 1:numel(Bt0) %find each tight section area
                Lt0_ = zeros(size(Lt0));
                Lt0_(Lt0==i)=1;
                Lt(:,:,s) = Lt0_;
                s = s + 1;
            end
        end
        
            
        waitbar(k/length(B),hw,'sectioning');
    end
    close(hw)
    hw = waitbar(0,'tight sectioning');
    for k = 1:length(Bt) % each tight section
        boundaryT = Bt{k};
        %boundary = boundaryT;
        genStructRec; %generates structures 'dr' to compare

        if 0
            stackWrite(Lt,'Lt.tiff')
            figure(9); imagesc(rb); axis image
            figure(2); imagesc(rbThick); axis image
            figure(3); imagesc(wbThick); axis image
            figure(4); imagesc(wb); axis image     
            figure(5); imagesc(tBoundsImg0); axis image   
            figure(7); imagesc(tBoundsImg); axis image
        end

        [Y,X]=findPos(R.*Lt(:,:,k));
        %figure(16);imagesc(Lt(:,:,k))
        %figure(16);imagesc(wb)
        %figure(16);imagesc(rb)
        dmin=[];
        szStruct(k) = sum(sum(Lt(:,:,k)));
        for i = 1:numel(Y)
           d_ = sqrt((boundaryT(:,1)-Y(i)).^2+(boundaryT(:,2)-X(i)).^2);
           dmin(i) = min(d_);
        end
        bw = 1;
        [Nr,~]=histcounts(drmin,'BinWidth',bw);
        skCoeff = 0.5; % rate of edge points
        figure(1);        
        if Nr(1)/sum(Nr)>=skCoeff % skinny structure
            plot(boundaryT(:,2), boundaryT(:,1), 'r', 'LineWidth', 2)
        elseif isempty(dmin) % empty false structure
            plot(boundaryT(:,2), boundaryT(:,1), 'g', 'LineWidth', 2)
        else
            plot(boundaryT(:,2), boundaryT(:,1), 'k', 'LineWidth', 2)
            
            d(ilast+1:ilast+numel(dmin)) = dmin;
            dr(ilastr+1:ilastr+numel(drmin)) = drmin;
            drA(ilastrA+1:ilastrA+numel(drAmin)) = drAmin;
            ilast=i+ilast;
            ilastr = i2+ilastr;
            ilastrA = i3+ilastrA;
        end
        waitbar(k/length(Bt),hw,'tight sectioning');
    end
    close(hw)
    hold off
    mag=1;
    set(gcf,'units','pixels','Position',[120,120,size(R,2)*mag,size(R,1)*mag]); 
    set(gca,'units','pixels','Position',[0,0,size(R,2)*mag,size(R,1)*mag]);
    imgFig = getframe(gcf);
    dataImg = imgFig.cdata; 
    imwrite(uint16(dataImg),fnameB);
    
    %% plot
    figure(2); % distance to the edge plot
    % 1/3 data
    bw=1; % binwidth
    [Nd,xed]=histcounts(d,'BinWidth',bw); % xe: edges
    xc = xed(1:end-1)+bw/2; % center pos
    plot(xc-bw/2,Nd);
    Ns = sum(Nd);
    xm = min(xed);
    
    
    hold on
    % 2/3 ref flat
    [Nr,xe]=histcounts(dr,'BinWidth',bw); % Nr :reference
    ixm = find(xe==xm); 
    Nsr = sum(Nr(ixm:end));
    Nr = Ns/Nsr*Nr;
    xcr = xe(1:end-1)+bw/2; % center pos
    plot(xcr-bw/2,Nr,'r')
    
    % 3/3 ref preBleach image
    [NrA,xe]=histcounts(drA,'BinWidth',bw); % Nr :reference
    ixm = find(xe==xm); 
    NsrA = sum(NrA(ixm:end));
    NrA = Ns/NsrA*NrA;
    xcrA = xe(1:end-1)+bw/2; % center pos
    plot(xcrA-bw/2,NrA,'k')
    hold off
    
    xlabel('distance to the edge')
    ylabel('#pxs')
    title('distribution of distance to the edge')
    grid
    legend('data','ref','image')
    set(gcf,'units','pixels','Position',[800,200,510,550]); 
    set(gca,'units','pixels','Position',[50,50,450,450]);
    imgFig = getframe(gcf);
    dataImg = imgFig.cdata; 
%imwrite(uint16(dataImg),fnameD);
    
    %% XLS
    NdXLS = nan(size(Nr));
    xcXLS = nan(size(Nr));
    NdXLS(1:numel(Nd)) = Nd;
    xcXLS(1:numel(Nd)) = xc;
    ix0 = find(xcr==xcrA(1));
    NraXLS = nan(size(Nr));
    xcraXLS = nan(size(Nr));
    NraXLS(ix0:ix0+numel(NrA)-1) = NrA;
    xcraXLS(ix0:ix0+numel(NrA)-1) = xcrA;    
    matARR = [xcXLS' NdXLS' xcr' Nr' xcraXLS' NraXLS'];
    save(fnameMAT,'matARR');
    
    disp(fname(4:ix_(1)-1));
    %pause
    
end

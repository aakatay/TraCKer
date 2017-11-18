function [Bt,Lt,rb] = findBoundaries(R,mnSz,t)
    rb=[];
    if t == 1 % tight boundaries
        iR = abs(im2bw(R,0)-1);

        iRc = conv2(iR,ones(2),'same');

        iRcs = iRc;
        iRcs(iRc < 4)=0;

        iRcsShft1 = zeros(size(iRcs));
        iRcsShft2 = iRcsShft1;
        iRcsShft3 = iRcsShft1;
        iRcsShft1(2:end,2:end) = iRcs(1:end-1,1:end-1);
        iRcsShft2(2:end,:) = iRcs(1:end-1,:);
        iRcsShft3(:,2:end) = iRcs(:,1:end-1);
        iRcsH = im2bw(iRcs + iRcsShft1+ iRcsShft2+ iRcsShft3,0);


        %% boundaries
        [Bt,Lt]=bwboundaries(abs(1-iRcsH),'noholes');

        rb = zeros(size(R));
        for k = 1:length(Bt) % each tight section
            boundary = Bt{k};

            %plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1)

            %find each section surface
            for b = 1:size(boundary,1) % draw the boundary on image
                rb(boundary(b,1),boundary(b,2))=1;
            end
        end

        if 0 
            figure(1)
            imagesc(iR); axis image
            figure(3)
            imagesc(iRcs); axis image
            figure(4)
            imagesc(abs(1-iRcsH)); axis image
            figure(5)
            imagesc(rb); axis image
        end        
    else
        %% FIND BOUNDARIES
        cnvSz = 5;
        convWin = ones(cnvSz);
        Rc = conv2(R,convWin,'same'); % convolved R
        [B,L]=bwboundaries(Rc,4,'noholes');
        %a=zeros(10); a(3:7,3:7)=1;a(3,3)=0; [b,l]=bwboundaries(a,4); plotBoundaries(b,l)

        %% FIND TIGHT BOUNDARIES    
        dr=[];Bt={};s=1;Lt=[];
        %if numel(B)>1, disp('multiple sections'); end;
        %hw = waitbar(0,'sectioning');
        for k = 1:length(B) % all sections
            boundary = B{k};

            %find each section surface
            rb = zeros(size(R));
            for b = 1:size(boundary,1) % draw the boundary on image
                rb(boundary(b,1),boundary(b,2))=1;
            end
            wb =double(watershed(rb,4)); wb=wb-1;wb(wb~=0)=1; 

            % find tighest boundary
            rbThick = conv2(rb,ones(3),'same');
            rbThick0 = rbThick;
            rbThick(rb>0)=0; 
            rbThick(wb==0)=0;    
            wbThick =double(watershed(rbThick,8));
            tBoundsImg = zeros(size(wbThick)); % tight boundary
            tBoundsImg(wbThick>=2)=1; tBoundsImg0 = tBoundsImg;
            tBoundsImg(rbThick>0)=0;
            [Bt0,Lt0]=bwboundaries(tBoundsImg,'noholes'); % tight boundaries
            % remove small structs
    %         sZsT = sum(reshape(Lt0,size(Lt0,1)*size(Lt0,2),size(Lt0,3)),1);
    %         ixRemove = sZsT<mnSz;
    %         Lt0(:,:,ixRemove)=[];
    %         Bt0 = Bt0(~ismember(1:size(Bt0,1), ixRemove)); 
    
            
            

            for i = 1:numel(Bt0) %find each tight section area
                Lt0_ = zeros(size(Lt0));
                Lt0_(Lt0==i)=1;
                if sum(Lt0_(:)) < mnSz, continue; end
                Lt(:,:,s) = Lt0_;
                Bt = [Bt;Bt0(i)];
                s = s + 1;
            end

            if 0 && ~isempty(Bt) % debug
                %plot(boundary(:,2), boundary(:,1), 'k', 'LineWidth', 2)
                %figure(11); imagesc(Lt0); axis image
                %figure(12); imagesc(Lt); axis image
    %            stackWrite(Lt,'Lt.tiff')
                saveFigImgFN = 'findBoundariesDBG.tif';
                figure(501); imagesc(im2bw(Rc,0)); axis image
                delete(saveFigImgFN);
                figure(1); imagesc(rb); axis image; title('rb (R\_5pxConv\_border)');saveFigImg;
                figure(2); imagesc(wb); axis image; title('wb (rb\_watershed)');saveFigImg;
                figure(3); imagesc(rbThick0); axis image; title('rbThick0 (rb\_3pxConv)');saveFigImg;
                figure(4); imagesc(rbThick); axis image; title('rbThick (rbThick0-rb-inv(wb))');saveFigImg;
                figure(5); imagesc(wbThick); axis image; title('wbThick (rbThick\_watershed)')   ;saveFigImg;
                figure(6); imagesc(tBoundsImg0); axis image; title('tBoundsImg0 (wbThick>=2)')  ;saveFigImg;
                figure(7); imagesc(tBoundsImg); axis image; title('tBoundsImg');saveFigImg;
                %figure(16);imagesc(Lt(:,:,k))
                figure(19); plotBoundaries(Bt,sum(Lt,3).*R,1); title('Bt.*R (R\_filtWith_tBoundsImg)');saveFigImg;
                ccc=3;
                
            end           
            %waitbar(k/length(B),hw,'sectioning');
        end
        %close(hw)
    end
    
    % remove single points
    ixDel = [];
    for k = 1:length(Bt)
        if size(Bt{k},1)==2 % single point
            ixDel = [ixDel k];
        end
    end
    ixDel = fliplr(ixDel);
    if ~isempty(ixDel)
        for k = 1:numel(ixDel)
            Bt{ixDel(k)}=[];
            Lt(:,:,ixDel(k))=[];
        end
        Bt= Bt(~cellfun('isempty',Bt));
    end
    
end







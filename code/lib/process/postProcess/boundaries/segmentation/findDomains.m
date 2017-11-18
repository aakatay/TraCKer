function [Rd, strInfo,Btb] = findDomains(Rc)
% find domains of a structure given recruitments
    %load Rc

    mbSz = 26; % min main body size
    sbSz = 7; % min side body size
    dbg = 0; % is debug 
    Btb = [];    
    if dbg, figure(88);imagesc(Rc); axis image; end

    %% BW rec

    RcBW = im2bw(Rc,0);
    if dbg, figure(89); imagesc(RcBW); axis image; end
    [szy,szx]=size(Rc);

    Rb = zeros(size(Rc));
    %%
    if 0
        figure(90);
        B = bwboundaries(RcBW,8,'noholes');
        axis image
        for i = 1:numel(B)
            if size(B{i},1)<=4
                continue;
            end

            c=cell2mat(B(i));
            plot(c(:,2),size(Rc,1)-c(:,1)+1,'LineWidth',1);
            axis image
            hold on
        end
        hold off
        xlim([1 szx])
        ylim([1 szy])
        axis image
    end
    %

    B = bwboundaries(RcBW,4);

    for i = 1:numel(B)
        if size(B{i},1)<=4
            continue;
        end

        c=cell2mat(B(i));
        if dbg, figure(91); plot(c(:,2),size(Rc,1)-c(:,1)+1,'LineWidth',1); axis image; hold on; end;
        for j= 1:size(c,1)
            Rb(c(j,1),c(j,2))=1; % boundary image
        end


    end

    if dbg, hold off; xlim([1 szx]);ylim([1 szy]); end;



    % boundary image
    if dbg, figure(92); imagesc(Rb); axis image; end;


    %% high res boundary image

    structHRfn = 'structHR.tif';
    RbH = zeros(size(Rc)*2);
    for i = 1:numel(B)
        if size(B{i},1)<=4
            continue;
        end

        c=cell2mat(B(i))*2;
        %plot(c(:,2),size(Rc,1)-c(:,1)+1,'LineWidth',1);
        c0x = 0;
        c0y = 0;
        for j= 1:size(c,1)
            dx = c0x-c(j,2);
            dy = c0y-c(j,1);
            if (abs(dx) + abs(dy)) == 2
                RbH(c(j,1)+dy/2,c(j,2)+dx/2)=1;
            end
            RbH(c(j,1),c(j,2))=1; % boundary image
            c0x = c(j,2);
            c0y = c(j,1);
            %imagesc(RbH)
            %pause
        end

    end

    if dbg, figure(93); imagesc(RbH);axis image; end;
    %% fill inside the boundaries

    RbHf = imfill(RbH);
    if dbg, figure(97); imagesc(RbHf); axis image; end;

    %% back to real resolution

    Rbf = RbHf([1:2:end]+1,[1:2:end]+1);
    if dbg, figure(98); imagesc(Rbf); axis image; end;

    %% struct numbers

    ws_ = watershed(RbH,8);
    if dbg, figure(99); imagesc(ws_); axis image; end;

    %% 
    ns = max(ws_(:))-1;
    if ~ns % one sparse domain
        Rd = double(im2bw(Rc,0));
        [Bt,Lt,rb] = findBoundaries(Rc,0,0);
        sz = sum(sum(sum(Lt)));
        strInfo = [sz 4 1]; % growth
        return
    end

    %% determine split structs and borders of struct
    cw = ones(3);
    for i = 1:ns
        L_ = zeros(size(RbH));
        L_(ws_==i+1)=1;
        L(:,:,i)=im2bw(conv2(L_,cw,'same'),0);
        szL(i) = sum(sum(L(:,:,i)));
    end
    [~,ixL]=sort(szL,'ascend');
    L = L(:,:,ixL);
    ws = ws_;
    ws(ws_==0)=1;
    ws(ws_==1)=0;
    ws=im2bw(ws,0);
    RbHs=im2bw(RbHf,0);


    sb = RbHs-ws; % structure branches
    if dbg, figure(100); imagesc(uint8(sb*50)+ws_); axis image; end;

    %% number branches

    sbH = conv2(repelem(sb,3,3),ones(3));
    sbH(sbH==9)=0;
    sbHws = watershed(sbH,8);
    sbws = double(sbHws(3:3:end,3:3:end)-1);
    nb = max(sbws(:));
    % colormap
    cmn=double(nb)+1;
    cmnr = 1:cmn;
    [~,ixs] = sort(rand(1,cmn));
    cix = cmnr(ixs);
    CM = jet(cmn);
    CM=CM(cix,:);
    CM(1,:)=0;
    if dbg, figure(1010101); colormap(CM); imagesc(sbHws); axis image; end;
    if dbg, figure(101); colormap(CM); imagesc(sbws); axis image; end;

        
    %% find connecting branches
    ps = [0 1 0; 1 0 1; 0 1 0]; % '+' conv window (plus sign)
    cnx = {};
    k=1;
    sbc2 = zeros(size(sbws));
    for i = 1:nb % each branch
        sbc = zeros(size(sbws));
        sbc(sbws==i) = 1; % current sb
        sbcArr(:,:,i) = sbc;
        sbcc = im2bw(conv2(sbc,ps,'same'),0); 
        cnx_ = [];
        for j = 1:ns % check each structure
            sbcx = sbcc+L(:,:,j);
            if ~isempty(find(sbcx==2)) % connected to a structure body
                cnx_ = [cnx_ j];
            end
        end
        cnx{i}={cnx_};
        cnxN(i) = numel(cnx_);
        if cnxN(i)>= 2
            cnxm{k} = cnx{i}; % multi connection
            sbc2 = im2bw(sbc,0)*k+sbc2; % structure branches (connecting)
            cbix(k)=i;
            k = k + 1;
            %imagesc(sbcRR+RcBW); axis image; insertBreakPoint
        end
    end
    
    
    %% join branches connected to connecting branches
    sbnc = im2bw(sbws,0) - im2bw(sbc2,0); % non-connecting branches
    sbnc(sbnc<0) = 0;
    sbnc = sbws.*sbnc;
    if dbg, figure(10000); imagesc(sbc2); axis image;end;
    %if dbg, figure(10000); imagesc(binImage(sbws,2)); axis image;end;
    
    %figure(10000); imagesc(sbc2); axis image
    %figure(101); colormap(CM); imagesc(sbws); axis image;
    %figure(134123); imagesc(Rc); axis image;
    
    ncb = max(sbc2(:));
    sbcRR = zeros(size(Rc));
    bccOv = zeros(size(Rc)); % overlap pos of two connecting branch
    sc = double(ws_) - 1; % structure cores
    sc(sc~=0) = 1;
    for i = 1:ncb % each connecting branch
        sbc_ = zeros(size(sbc2));
        sbc_(sbc2==i) = 1;
        sbcc = conv2(sbc_,ones(3),'same');
        sbcc = im2bw(sbcc,0);
        for j = 1:nb % find connected branches
            if cbix(i) == j, continue; end;
            sbnc2 = zeros(size(sbws)); % current sb
            sbnc2(sbws==j) = 1;     
            sbcOv = sbcc+sbnc2; % overlap
            if ~isempty(find(sbcOv==2)) % connected in HR
                sbc_rr = double(im2bw(binImage(sbc_,2),0));
                sbc_rrc = conv2(sbc_rr,ps,'same');
                sbc_rrc = im2bw(sbc_rrc,0);
                sbnc2rr = im2bw(binImage(sbnc2,2),0);
                sbcOv2 = sbc_rrc+sbnc2rr; % overlap
                if ~isempty(find(sbcOv2==2)) % connected in RR
                    if ~isempty(find(ismember(cbix,j)~=0)) % connected branch is a connecting branch
                        b1 = double(sbws==cbix(i));
                        b2 = double(sbws==j);
                        b1c = double(im2bw(conv2(b1,ps,'same'),0));
                        b2c = double(im2bw(conv2(b2,ps,'same'),0));
                        b12 = double((b1c + b2c)==2);
                        if dbg, figure(12121); imagesc(b12); axis image;end;
                        bccOv = bccOv + binImage(b12.*sc,2);
                    else
                        sbc_ = sbnc2+sbc_;
                    end                    
                end
            end
        end
        sbcRR = sbcRR + im2bw(binImage(sbc_,2),0)*i; % real resolution
    end
    sbcRR(bccOv>0)=0; 

    if dbg, figure(10101); imagesc(sbcRR); axis image;end;
    %% rec without branches 
    if ~exist('sbcRR'), 
        RcBWwob =  RcBW;
        ncb = 0; % # of connecting branches
    else
        RcBWwob =  RcBW - im2bw(sbcRR,0); % without branches
        ncb = max(sbcRR(:)); % # of connecting branches
    end
    
    if dbg, figure(102);imagesc(RcBWwob); axis image;end;
    %imagesc(sum(sbcRR,3))


    %% structure boundaries
    mnBsz = 4;
    mnBsz = 2;
    B = bwboundaries(RcBWwob,4,'noholes');
    Rb0 = zeros(size(Rc));
    if dbg,figure(191); axis image; end;
    for i = 1:numel(B)
        if size(B{i},1)<=mnBsz
            continue;
        end

        c=cell2mat(B(i));
        if dbg, plot(c(:,2),size(Rc,1)-c(:,1)+1,'LineWidth',1);axis image;hold on;end;
        for j= 1:size(c,1)
            Rb0(c(j,1),c(j,2))=1; % boundary image
        end

    end
    if dbg, hold off; xlim([1 szx]);ylim([1 szy]); end;


    % boundary image
    if dbg, figure(192); imagesc(Rb0);axis image;end;
    %% high res boundary image
    RbH0 = RbH;
    RbH = zeros(size(Rc)*2);
    for i = 1:numel(B)
        if size(B{i},1)<=mnBsz
            continue;
        end

        c=cell2mat(B(i))*2;
        %plot(c(:,2),size(Rc,1)-c(:,1)+1,'LineWidth',1);
        c0x = 0;
        c0y = 0;
        for j= 1:size(c,1)
            dx = c0x-c(j,2);
            dy = c0y-c(j,1);
            if (abs(dx) + abs(dy)) == 2
                RbH(c(j,1)+dy/2,c(j,2)+dx/2)=1;
            end
            RbH(c(j,1),c(j,2))=1; % boundary image
            c0x = c(j,2);
            c0y = c(j,1);
            %imagesc(RbH)
            %pause
        end
    end
    if dbg, figure(193); imagesc((RbH));axis image;end;
    if dbg, figure(193193); imagesc((RbH0));axis image; end;
    %% fill inside the boundaries

    RbHf = imfill(RbH);
    RbH0f = imfill(RbH0);
    if dbg, figure(197);imagesc(RbHf); axis image; end;
    if dbg, figure(197197);imagesc(RbH0f); axis image; end;
    RbHfrr = RbHf([1:2:end]+1,[1:2:end]+1);
    RbH0frr = RbH0f([1:2:end]+1,[1:2:end]+1);
    RbMiss=(RbH0frr-RbHfrr-im2bw(sbcRR,0)).*Rc>0;
    if dbg, figure(197197197);imagesc(RbMiss); axis image; end;
    %% back to real resolution
    Rbf = RbHf([1:2:end]+1,[1:2:end]+1) + RbMiss*2;
    if dbg, figure(198); imagesc(Rbf); axis image; end

    %% fill each struct with different color
    RbCfull = zeros(size(Rb0));
    RbCarr = zeros(size(Rb0));
    ns = numel(B);
    k=1;
    strInfo=[];
    for i = 1:ns
        Rb = zeros(size(Rb0));
        if size(B{i},1)<=mnBsz
            continue;
        end
        c=cell2mat(B(i));
        for j= 1:size(c,1)
            Rb(c(j,1),c(j,2))=1; % boundary image
        end
        %imagesc(Rb); pause
        RbC = imfill(Rb);
        RbCfull = RbC*k+RbCfull;

        % label structures
        strInfo(k,1) = sum(RbC(:)>0); % structure size

        k = k +1;
        %imagesc(RbCfull); axis image; pause
        if ~isempty(find(RbCfull>i))
            imagesc(RbCfull); axis image
            pause
        end
    end
    %% add missing structures
    RbMissCv = conv2(double(RbMiss),ones(2),'same'); 
    RbMissPeak = imregionalmax(RbMissCv);
    nMiss = numel(find(RbMissPeak>0));
    RbMissPeak(RbMissPeak>0)=1:nMiss;
    if dbg, figure(197197198); imagesc(RbMissPeak); axis image; end
    for i = 1: nMiss % each missing structure
        RbMissPeak_ = zeros(size(RbMissPeak));
        RbMissPeak_(RbMissPeak == i) = 1;
        RbMissPeakCv = conv2(RbMissPeak_,ones(3),'same');
        RbMiss_ = RbMiss.*RbMissPeakCv;
        RbCfull = RbCfull + RbMiss_*(max(RbCfull(:))+i);
        strInfo(end+1,1) = sum(RbMiss_(:));
    end
    %%
    cmn=max(RbCfull(:))+1;
    CM = jet(cmn);
    cmnr = 1:cmn;
    [~,ixs] = sort(rand(1,cmn));
    cix = cmnr(ixs);
    CM=CM(cix,:);
    CM(1,:)=0;
    if dbg, figure(201); colormap(CM); imagesc(RbCfull); axis image; end


    % label structures
    %strInfo(:,2) = 11*(strInfo(:,1)>=mbSz); % (side) bodies (lb:11)
    strInfo(:,2) = 11; % (side) bodies (lb:11)
    [~,ixmx] = max(strInfo(:,1));
    strInfo(ixmx,2)  = 1; % main body (lb:1)
    strInfo(:,3) = 1:size(strInfo,1); % initialie group indices



    %% add connecting branches
    if dbg,figure(202); imagesc(sbcRR); axis image; end;

    RbCfullwB = RbCfull;
    if dbg, figure(1110);imagesc(RbCfull); axis image; end; % with branches
    bpn=[];
    msix = find(strInfo(:,2)>0); % main structure index
    nms = numel(msix); % # of main structs
    ns = max(RbCfull(:)); % # of structs
    cnix = cell(ncb,1); % indices of structures connected with the cb
    for i = 1:nms % each main body
        if ~ncb, break; end; % no branches
        msix2 = msix(i);
        ms = double(RbCfull==msix2); % main structure
        ixExc = msix2;
        while 1 % loop till no connection is found
            msc = im2bw(conv2(ms,ps,'same'),0); 
            cx = 0;
            for j = 1:ncb % each connecting branch
                cb = zeros(size(sbcRR));
                cb(sbcRR==j) = 1; % connecting branch
                if ~isempty(find((cb + msc)>=2)) % connected to main struct

                    for k = 1:ns % find which side structs are connected
                        if ismember(k,ixExc),
                            continue; 
                        end; % skip if main struct or connected side struct
                        cbc = im2bw(conv2(cb,ps,'same'),0); 
                        ss= double(RbCfull==k); % side structure
                        if ~isempty(find((cbc + ss)>=2)) % connected to side struct
                            strInfo(k,2)=round(strInfo(msix2,2)/10)*10+2; % connected structs (rem(lb,10)=2)
                            strInfo(k,3)=strInfo(msix2,3); % group assignment
                            ms = ms + ss;
                            cx = 1; % connection found
                            RbCfullCurr = RbCfull;
                            RbCfullCurr(~ismember(RbCfull,ixExc))=0;
                            RbCfullCurrBW = im2bw(RbCfullCurr,0);
                            RbCfullCurrOv = RbCfullCurrBW+cbc;
                            msixCurr=unique(RbCfullCurr(find(RbCfullCurrOv==2)));
                            cnix{j} = [msixCurr, k];
                            ixExc = [ixExc' k]';
                        end % #if connected to the connecting branch
                    end % #for each side struct
                end % #if connected to main struct
            end % #for each connecting branch
            if cx == 0, break; end; % no connection is found
        end % #while till no connection
        %RbCfullwB = RbCfullwB + sbcRR(:,:,i)*(1+max(RbCfullwB(:)));%double(cnxs(kk));
        %insertBreakPoint
    end
    if dbg, figure(1111);imagesc(RbCfullwB); axis image; end; % with branches
    %figure(1001); imagesc(ismember(RbCfullwB,[11 4])); axis image
    
    % correct main body strType#
    [~,ixmx]=max(strInfo(:,1));
    ngm = find(strInfo(:,3)==strInfo(ixmx,3)); % group of the main structure
    ngm(ngm==ixmx)=[]; % other than main struct
    strInfo(ixmx,2) = 1;
    strInfo(ngm,2) = 2;
    %% assign branch points to structures
    for i = 1:ncb % each connecting branch
        bpn = [];
        [y,x] = find(sbcRR==i);
        for j = 1:numel(y) % each point on the branch
            bp = zeros(size(sbcRR(:,:,1)));
            bp(y(j),x(j))=1;
            bpc = im2bw(conv2(bp,ps,'same'),0);
            cnxs = cnix{i};
            for k = 1:numel(cnxs) % each connected structure to the branch
                s = double(RbCfullwB==cnxs(k));
                bpcx = double(bpc)+s;
                bpn(j,k) = numel(find(bpcx>=2)); % #neighb of bp
            end
        end
        bpn0 = bpn;

        for j = 1:numel(cnxs)
            cnxsSz(j) = numel(find(RbCfullwB==cnxs(j)));
        end
        if min(cnxsSz)<sbSz % too small structs
            [~,ixmn] = min(cnxsSz);
            [~,ixmx] = max(cnxsSz);
            ixmn = cnxs(ixmn);
            ixmx = cnxs(ixmx);
            RbCfullwB(RbCfullwB==ixmn) = ixmx;
            RbCfullwB(sub2ind(size(RbCfullwB),y,x)) = ixmx;
            RbCfullwB(RbCfullwB>ixmn)=RbCfullwB(RbCfullwB>ixmn)-1;
            strInfo(ixmn:end-1,:)=strInfo(ixmn+1:end,:);
            strInfo(end,:)=[];
            strInfo(strInfo(:,3)>ixmn,3)=strInfo(strInfo(:,3)>ixmn,3)-1;
            for j =1:numel(cnix)
                cnix{j}(cnix{j}==ixmn) = ixmx;
                cnix{j} = cnix{j}-(cnix{j}>ixmn);
            end
            continue;
        end
        [py,sx]=find(bpn>=2); % point ix(y) and structure ix(x)
        [py, pix ] = sort(py);
        sx = sx(pix);

        RbCfullwBnp = zeros(size(RbCfullwB));
        while ~isempty(py) % assign 2-edge connected points
            RbCfullwB(y(py(1)),x(py(1))) = cnxs(sx(1));
            RbCfullwBnp(y(py(1)),x(py(1))) = cnxs(sx(1)); % new points
            bpn(py(1),:) = [];
            y(py(1)) = []; x(py(1)) = [];
            py(1)=[]; sx(1)=[];
            py = py - 1;
        end
        [~,mx]=max(cnxsSz); % largest structure
        isCx = 0;
        while 1 % lookfor points connected to largest structure
            isCx = 0;
            ixDel = [];
            for j = 1:size(bpn,1) % each point on the branch
                sp = zeros(size(RbCarr(:,:,1)));
                RbCfullwBmx = zeros(size(RbCarr(:,:,1)));
                RbCfullwBmx(RbCfullwB == cnxs(mx)) = 1; % largest structure
                sp(y(j),x(j)) = 1; % single point
                spc = conv2(sp,ps,'same'); % conv with '+'
                if ~isempty(find(RbCfullwBmx+spc>=2)) % connected
                    RbCfullwB(y(j),x(j)) = cnxs(mx);
                    isCx = 1;
                    ixDel = [ixDel j];
                end
            end
            bpn(ixDel,:) = []; y(ixDel) = []; x(ixDel) = [];
            if ~isCx
                break; 
            end
        end
        while ~isempty(bpn) % if more points on the branch
            % check for connection
            for j = 1:size(bpn,1)
                sp = zeros(size(RbCarr(:,:,1)));
                sp(y(j),x(j)) = 1; % single point
                spc = conv2(sp,ps,'same'); % conv with '+'
                spOvLap = RbCfullwB + spc; % overlap with structures
                spOvLapBW = im2bw(RbCfullwB,0) + spc; % overlap with structures
                [ys,xs] = find( spOvLapBW >=2); %overlapping structure points
                if ~isempty(ys)
                    RbCfullwB(y(j),x(j)) = spOvLap(ys,xs)-1; % is of overlapping struct
                    bpn(j,:) = []; y(j) = []; x(j) = [];
                end
            end
        end    
    end
    if dbg, figure(1112);imagesc(RbCfullwB); axis image; end; % with branches
    %% clean strInfo
    %strInfo(find(strInfo(:,2)==0),:)=[];

    %% display groups in same color
    gixs = strInfo(:,3); % group indices
    gix = unique(gixs,'first');
    RbCfullwBgr = RbCfullwB;
    for i = 1:numel(gix)
        [ixs] = find(strInfo(:,3) == gix(i));  % structure ix
        RbCfullwBgr(ismember(RbCfullwBgr,ixs))=gix(i);
        %imagesc(ismember(RbCfullwBgr,ixs))
    end
    if dbg,figure(1113);imagesc(RbCfullwBgr); axis image; colormap('jet');end; % with branches
    ixmb=find(ismember(strInfo(:,2),[1,2,11,12]));
    RbCfullwBmb = zeros(size(RbCfullwB)); % main bodies
    RbCfullwBmb(ismember(RbCfullwBgr,ixmb))=1;
    if dbg,figure(1114);imagesc(RbCfullwBmb); axis image;end; % only bodies
    strInfo0=strInfo;

    %% remove other than bodies (strInfo & RbCfullwB)
    strInfo = strInfo0;
    ixDel=find(strInfo(:,2)==0);
    ixKeep=strInfo(find(strInfo(:,2)~=0),3);
    ix = 1:size(strInfo,1);
    ixDel0 = ixDel;

    for j=1:numel(ixDel)
        i=ixDel(j);
        ixKeep(ixKeep>i)= ixKeep(ixKeep>i)-1;
        ixDel(j+1:end)=ixDel(j+1:end)-1;
    end
    strInfo(strInfo(:,2)==0,:)=[];
    strInfo(:,3)=ixKeep;
    ixDel = ixDel0;
    while ~isempty(ixDel)
        i = ixDel(1);
        RbCfullwB(RbCfullwB==i) = 0;
        RbCfullwB(RbCfullwB>i)= RbCfullwB(RbCfullwB>i)-1;
        ixDel(1)=[];
        ixDel = ixDel - 1;
    end
    if dbg,imagesc(RbCfullwB); axis image; end

    
    %% update small bodies
    [ixBody,~] = find(ismember(strInfo(:,2),[1,2,11,12]));
    for i = 1:numel(ixBody)
        RcBody = zeros(size(Rc));
        RcBody(RbCfullwB==ixBody(i))=1;
        [Bt,Lt,rb] = findBoundaries(RcBody,0,1);
        sz=sum(Lt(:)>0);
        if sz < sbSz % too small
            ng = find(strInfo(:,3)==strInfo(ixBody(i),3));
            if numel(ng)>1
                strInfo(ng,3)=strInfo(ng(2),3);
            end
            strInfo(ixBody(i),:)=[];    
            strInfo(strInfo(:,3)>ixBody(i),3)=strInfo(strInfo(:,3)>ixBody(i),3)-1;
            ixBody(i+1:end)=ixBody(i+1:end)-1;
            RbCfullwB(RbCfullwB==ixBody(i))=0;
            RbCfullwB(RbCfullwB>ixBody(i))=RbCfullwB(RbCfullwB>ixBody(i))-1;
        else
            strInfo(ixBody(i),1)=sz;
        end
        ccc=3;
    end
    
    RbCfullwB0=RbCfullwB;
    %% add edge recruits
    RbCfullwB=RbCfullwB0;

    mnSz=0;
    cvSz = 5;
    Ng = unique(strInfo(:,3)); % # of bodies
    nb = numel(Ng);
    ns = max(RbCfullwB(:));
    RcEdgeAll = zeros(size(RbCfullwB));
    for i = 1:nb % each main body
        ng = find(strInfo(:,3)==Ng(i));
        Rg = RbCfullwB; 
        Rg(find(ismember(RbCfullwB,ng)==0))=0; % selected group recs
        RgCv = im2bw(conv2(Rg,ones(cvSz),'same'),0);
        RgCvB = RgCv - im2bw(RbCfullwB,0);
        RgCvB(RgCvB<0)=0;
        RcEdge = im2bw(Rc.*RgCvB,0);
        [y,x] = find(RcEdge>0);
        for j = 1:numel(y) % if connected to a body
            RcEdgeSp = zeros(size(RcEdge));
            RcEdgeSp(y(j),x(j))=1;
            RcEdgeSpCv = conv2(RcEdgeSp,ps,'same');
            sumBW = RcEdgeSpCv + im2bw(RbCfullwB,0);
            if ~isempty(find(sumBW==2)) % add to the body
                ixBd = unique(RbCfullwB(sumBW==2));
                if numel(ixBd)>1, error('needs further coding'); end
%RbCfullwB(y(j),x(j)) = ixBd;
%RcEdge(y(j),x(j))=0;
            end
        end
        RcEdgeAll = RcEdgeAll+RcEdge;
        strTyp = 10*round(min(strInfo(ng,2))/10)+3; % structure type : edge (rem(lb,10)=3)
        sz = sum(RcEdge(:));
        if ~sz, 
            continue; 
        end;
        ns = ns + 1;
        strInfo(ns,:)=[sz strTyp i];
        RbCfullwB = RbCfullwB + RcEdge*ns;
    end

    if dbg, figure(1115);imagesc(RbCfullwB); axis image; colormap('jet'); end; 
    %% add grow recruits
    RcEdge = RcEdgeAll;
    RcGrow00 = Rc-Rc.*im2bw(RbCfullwB,0);
    if dbg, figure(1116); imagesc(RcGrow00); axis image; end;
    cvSz = 7;
    RcGrow0 = RcEdge.*im2bw(conv2(RcGrow00,ones(cvSz),'same'),0)+RcGrow00; % add recs from edge
    RcGrow = RcEdge.*im2bw(conv2(RcGrow0,ones(3),'same'),0)+RcGrow0; % add immediate neighbours
    if dbg, figure(1117);imagesc(RcGrow); axis image; end

    RcMain = logical(RbCfullwB)-logical(RcGrow);
    RcMain(RcMain<0)=0;
    [Btb,Ltb,rb] = findBoundaries(RcGrow,mnSz,0);

    if dbg, figure(1118); plotBoundaries(Btb,RcMain+sum(Ltb,3).*Rc); end;
    RcDom = RbCfullwB.*RcMain.*im2bw(Rc,0); % domains
    ns = size(strInfo,1);
    RcDomSp=zeros(size(RcDom));
    if ~isempty(Ltb)
        for i = 1:size(Ltb,3)
            sz = sum(sum(Ltb(:,:,i)));
            if sz == 1, 
                RcDomSp = RcDomSp + Ltb(:,:,i).*RbCfullwB.*im2bw(RcGrow,0);
                continue; 
            end;
            ns = ns + 1;
            strInfo(ns,:) = [sz 4 0]; % (lb:4)
            RcDom = RcDom + im2bw(Ltb(:,:,i).*RcGrow,0)*ns;
        end
    end
    RcDom = RcDom + RcDomSp;
    if dbg, figure(1119); imagesc(RcDom); axis image; colormap('jet'); end

    %if dbg, figure(1120); imagesc(Rc); axis image; end

    %% update strInfo and RcDom
    
    % update edge sizes
    ixEdge = find(ismember(strInfo(:,2),[3,13]));
    for i = 1:numel(ixEdge)
        sz = sum(sum(RcDom==ixEdge(i)));
        strInfo(ixEdge(i),1)=sz;
    end

    %% add single points
    sp = RcBW-im2bw(RcDom,0); % single points
    sp(sp<0)=0;
    Rd = RcDom; % domains
    nRd = double(max(Rd(:)))+1;
    Rd = Rd + sp.*(nRd);
    
    nRd = nRd+1;
    if dbg,
        figure(1121);
        imagesc(Rd); 
        axis image; 
        % colormap
        cmn = nRd+1;
        cmnr = 1:cmn;
        CM = jet(cmn);
        [~,ixs] = sort(rand(1,cmn));
        cix = cmnr(ixs);
        CM=CM(cix,:);
        CM(1,:)=0;
        colormap(CM); 
    end; % with branches
    
    
    
    
    
    
    
    if ~dbg, return; end;

    figure(111212121); imagesc(Rc); axis image;
    figure(1121);
    pause;
    
    return
    % edges
    if dbg, figure(1122); imagesc(Rd.*ismember( Rd, find(ismember(strInfo(:,2),[3,13])) )); axis image; end
    % growth
    if dbg, figure(1123); imagesc(Rd.*ismember( Rd, find(ismember(strInfo(:,2),[4])) )); axis image; end
    % main body
    if dbg, figure(1124); imagesc(Rd.*ismember( Rd, find(ismember(strInfo(:,2),[1,11])) )); axis image; end
    % side body
    if dbg, figure(1125); imagesc(Rd.*ismember( Rd, find(ismember(strInfo(:,2),[2,12])) )); axis image; end
    % single points
    if dbg, figure(1126); imagesc(Rd==max(Rd(:))); axis image; end
    %% plot strInfo

    figure(1130);
    for i = 1:size(strInfo,1)
        RdD = Rd==i;
        imagesc(RdD); axis image; 

        disp(sprintf('strInfo : %3i \t %3i \t %3i',strInfo(i,:)));
        pause 
    end
    %%
    figure(1131);
    ngs = unique(strInfo(:,3));
    nG = numel(ngs);
    RdDs = zeros(size(Rc)); % sum
    for i = 1:nG
        ng = find(strInfo(:,3)==ngs(i));
        RdD = ismember(Rd,ng);
        imagesc(RdD); axis image; 
        RdDs = RdDs +RdD;

        %disp(sprintf('strInfo : %i \t %i \t %i',strInfo(i,:)));
        pause

    end
    if dbg, figure(1132); imagesc(RdDs); axis image; end

    if dbg, figure(1120); imagesc(Rc); axis image; end


    %%

    if 0 
    %% sort figures
        nf=[];
        figHandles = findobj('Type','figure');
        for i = 1:numel(figHandles)
            nf(i) = figHandles(i).Number;
        end
        nf=sort(nf,'ascend');
        nf=sort(nf,'descend');
        for i = 1:numel(figHandles)
            figure(nf(i));
        end
    %%
    end

end
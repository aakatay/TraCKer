%% 1- select a rectangle to remove outside dots (stay in image size)
%% 2- to quit select a rect towards outside on upper left of the image 
%% 3- it will find edge and distance to the edge with reference and plot and write to xls
clear
%close all
isAutoBorderDetection = 1;
if isempty(strfind(pwd,'recComp'))
    mkdir('recComp')
    cd('recComp')
end
cropTimeWin; % check drift, find a time window(f1,f2)
fPreImg.name = sprintf('..\\%s',fimg.name);

isAutoFindDomains = 0;
selMask = -1; % negative mask
selMask = 1;
selMask = 0; % no cropping
mnSz0 = 65; % only large structures are processed
%mnSz0 = 20; % only large structures are processed
mnDensity = 0.25; % only dense recruiting structs are processed

mnSz = 0; % only large structures are processed
fnameRecSumDIR = rdir('..\binImgRcrtSum_time*');
fnameRecTimeDIR = rdir('..\binImgRcrtTime_time*');
fnameRecSumBlchDIR = rdir('..\binImgRcrtSum_stBlchTime*');
fnameRecTimeBlchDIR = rdir('..\binImgRcrtTime_stBlchTime*');
fnameLapMovie = rdir('..\lap_*'); % pre bleach image
fnameMask = 'posMask.mat';
cfg = struct;
cfg.process = struct;
edgeDistDataDIR = rdir('edgeDistData*.mat');
isSelSt = 0;
if ~isempty(edgeDistDataDIR), isSelSt=1; end; % edge analysis done display specs (different mode)
% file label
dirName = pwd; % directory name
dS = strfind(dirName,'\'); % dirSlash
cellName = dirName(dS(end-1)+1:dS(end)-1);
dateName = dirName(dS(end-2)+1:dS(end-2)+6);
if selMask
    if selMask == 1
        maskLab = '-mask1'; % mask label
    elseif selMask == -1
        maskLab = '-mask2'; % mask label
    else
        error('seleck mask')
    end
else
    maskLab = ''; 
end
    
cellLabel = [dateName '-' cellName maskLab];

%% structure replacement
    %[intCoef4Lap,smCoeff] = findDiffRec;

%% laptime dynamics
    
    %noiseThresh = smCoeff*intCoef4Lap;
    noiseThresh =0;
    
    isNoLapAcq = 0;
    if isempty(fnameLapMovie)
        fnameLapMovie = fPreImg;
        warning('no pre-lap-acquisition');
        isNoLapAcq = 1;
    end
    [Adev,Afamp,Amean] = findFlickering(fnameLapMovie.name,noiseThresh);
    A = repelem(Amean,4,4);
    Afamp = repelem(Afamp,4,4);
    Adev = repelem(Adev,4,4);

%%
ps= []; % processed structs types (for color code)

close all;
fnameRecSum = fnameRecSumDIR(1).name;
fnameRecTime = fnameRecTimeDIR(1).name;
stMapImgFN= [ 'structMap' cellLabel '.tif']; % map
stMapUpdImgFN= [ 'structMapUpd' cellLabel '.tif']; % map
stMapUpdImgFN2= [ 'structMapUpdSelSt' cellLabel '.tif']; % map
stMapMATFN= [ 'structMap' cellLabel '.mat']; % map
stMapHRMATFN= [ 'structMapDomainsHR.mat']; % map
stMapUpdMATFN= [ 'structMapUpd' cellLabel '.mat']; % map
overlayFullFN{1} = [ 'overlayFull_' cellLabel '.tif']; % map
overlayFullFN{2} = [ 'overlayFull2_' cellLabel '.tif']; % map
overlayFullFN{3} = [ 'overlayFull3_' cellLabel '.tif']; % map
R = double(imread(fnameRecSum));
Rtime = double(imread(fnameRecTime));

%% time normalization
Clevel = 2^6;
Rtime = Rtime*(Clevel-1)/max(Rtime(:));

%% high resolution overlay image
dc = 1; % data crop #
mag = 2;
[Y,X] = find(R>0);
Y = Y*mag-floor(mag/2);
X= X*mag-floor(mag/2);
XY = (X-1)*size(R,1)*mag+Y;
O = zeros(size(R));
Clevel = 2^16;
CM = parula(Clevel);

if ~isNoLapAcq
    A1 = repelem(A,mag,mag);
    A2 = repelem(Afamp,mag,mag);
    A3 = repelem(Adev,mag,mag);

    A123 = cat(3,A1,A2);
    A123 = cat(3,A123,A3);

    for i = 1:numel(overlayFullFN)
        Ap = A123(:,:,i);
        Amag = double(Ap);
        Amag = round(Amag*(Clevel-1)/max(Amag(:)))+1;
        Acol1D = CM(Amag(:),:); 
        Acol1D = reshape(Acol1D,size(Amag,1),size(Amag,2),3);

        % red points
        Acol1Dch1 = Acol1D(:,:,1); % red 
        Acol1Dch2 = Acol1D(:,:,2);
        Acol1Dch3 = Acol1D(:,:,3); % blue
        Acol1Dch1(XY)=1;
        Acol1Dch2(XY)=0;
        Acol1Dch3(XY)=0;
        Acol1D(:,:,1) = Acol1Dch1;
        Acol1D(:,:,2) = Acol1Dch2;
        Acol1D(:,:,3) = Acol1Dch3;

        AcolImg = uint16((2^16-1)*Acol1D);
        imwrite(AcolImg,overlayFullFN{i},'Compression','none')
    end
end

%%  define a mask to process in multiple steps
if ~selMask
    ; % no mask selection
elseif ~exist(fnameMask)
    mag = 2;
    figure(41)
    imagesc(A)
    hold on
    [Y,X] = find(R>0);
    scatter(X,Y,'.','r')
    hold off
    set(gcf,'units','pixels','Position',[720,-120,size(R,2)*mag,size(R,1)*mag]); 
    set(gca,'units','pixels','Position',[0,0,size(R,2)*mag,size(R,1)*mag]);
    % draw a poly for separating the image two pieces
    hp = impoly(gca);
    pos = getPosition(hp);

    % move frist and last points to the edge
    pos(1,2) = 0;
    pos(end,2) = size(A,1);

    pos1 = zeros(size(pos,1)+2,2); % crop 1
    pos1(2:end-1,:) = pos;
    pos1(1,:) = [0,0];
    pos1(end,:) = [0,size(A,1)];
    mask = roipoly(A,pos1(:,1),pos1(:,2));
    imagesc(mask)
    save(fnameMask,'pos1','mask');
else 
    load(fnameMask);
end

if selMask
    A = A.*uint16(mask);
    R = R.*double(mask);
end


%% =========== find boundaries =============
mag=1;
if isCropTimeWin | (~exist(stMapMATFN)  & ~exist(stMapUpdMATFN)) % find coarse boundaries
    %R = double(im2bw(R,0));
    mnFilt = [mnSz0 mnDensity];
    [Bt,Lt] = findBoundaries(R,mnFilt,0);
    
    tightPos = cell(size(Lt,3),1); % selected tight requests
    save(stMapMATFN,'Lt','Bt','tightPos','-v7.3') % structMap

    dispBoundaries(R,Bt,Lt,mag,[],isSelSt);
    imgFig = getframe(gcf);
    dataImg = imgFig.cdata; 
    imwrite(uint16(dataImg),stMapImgFN); % structMap
elseif exist(stMapUpdMATFN) % load updated file
    load(stMapUpdMATFN);
    Lt = LtNew;
    Bt = BtNew; clear BtNew LtNew;
    dispBoundaries(R,Bt,Lt,mag,ps,isSelSt);
        
else % exist(stMapMATFN)
    load(stMapMATFN);
    dispBoundaries(R,Bt,Lt,mag,ps,isSelSt);
end % already mapped once        

%% ===== process & display ====================================    
cfg.process.auto = 0;
if isAutoFindDomains % automatic clustering
    if exist(stMapHRMATFN)% display
        RdMb = zeros(size(R));
        RdcMb = zeros(size(R));
        load(stMapHRMATFN); % (Rdcrop strInfo BtbArr BtbIx pos)
        ns = numel(Rdcrop);
        for i = 1:ns
            Rdc = Rdcrop{i};
            stInf = strInfo{i};
            p = pos(i,:);
            %Btb = BtbArr{i};
            RdcMb(p(1):p(2),p(3):p(4)) = im2bw(Rdc,0)*i; % structure index
            RdMb(p(1):p(2),p(3):p(4)) = ceil(Rdc/max(Rdc(:))*64);
        end
        figure(11);
        plotBoundaries(BtbArr,RdMb,1);
        mag = 1;
        set(gcf,'units','pixels','Position',[120,120,size(R,2)*mag,size(R,1)*mag]); 
        set(gca,'units','pixels','Position',[0,0,size(R,2)*mag,size(R,1)*mag]);
        %imagesc(RdMb); axis image;

    else % call findDomains
        cfg.process.auto = 1;
        ns = size(Lt,3);
        Rd = zeros(size(R));
        RdMb = zeros(size(R));
        BtbArr = [];
        Bt = [];
        Lts = Lt;
        Lt = zeros(size(R));
        hw = waitbar(0);
        k=1;
        for i = 1:ns
            %if i == 34, continue; warning('skipping # 34'); end;
            %i = 34
            %i=4
            %i = 19;

            Rc = Lts(:,:,i).*R;
            [y_,x_]=find(Rc>0);   
            [ysz,xsz]=size(Rc);
            fs = 5; % frame size
            y1=min(y_)-fs;
            y2=max(y_)+fs;
            x1=min(x_)-fs;
            x2=max(x_)+fs;
            if y1<0, y1=1; end;  if y2>ysz, y2=ysz; end;
            if x1<0, x1=1; end;  if x2>ysz, x2=xsz; end;
            pos(i,:) = [y1 y2 x1 x2];
            Rs = Rc(y1:y2,x1:x2);
            [Rdcrop{i},strInfo{i},BtbCrop] = findDomains(Rs); % ====================== domain
            for j = 1:numel(BtbCrop)
                BtbC = BtbCrop{j};
                BtbC = [BtbC(:,1)+y1-1 BtbC(:,2)+x1-1];
                BtbCrop{j} = BtbC;
            end
            Rdc = Rdcrop{i};
            stInf = strInfo{i};
            Rd_(y1:y2,x1:x2) = ceil(Rdc/max(Rdc(:))*64);
            BtbArr = [BtbArr; BtbCrop];
            BtbIx{i} = numel(BtbArr)-numel(BtbCrop)+1:numel(BtbArr);
            % main body only
            ixmb = find(ismember(stInf(:,2),[1,2,11,12])); % main body indices
            if isempty(ixmb), continue; end;
            RdcMb = Rdc;
            RdcMb(~ismember(Rdc,ixmb))=0;
            RdMb(y1:y2,x1:x2) = ceil(RdcMb/max(RdcMb(:))*64);
            RdMb_ = zeros(size(RdMb));
            RdMb_(y1:y2,x1:x2) = RdMb(y1:y2,x1:x2);
            if 0  % disp the structure
                figure(11);
                imagesc(Rdcrop{i}); axis image;
                if ~isempty(BtbCrop)
                    hold on 
                    plotBoundaries(BtbCrop,Rdc); 
                    hold off
                end
                figure(12);
                imagesc(A(y1:y2,x1:x2)); axis image;
                pause
            end
            cc=3;
            %% main body boundaries
            isTight = 1;
            nscixs = unique(RdMb_);
            nscixs(nscixs==0) = [];
            for j=1:numel(nscixs)
                nscix = nscixs(j);
                RdMb__ = RdMb_;
                RdMb__(RdMb_~=nscix)=0;
                mnSz = 15;
                [Btc,Ltc,rb] = findBoundaries(RdMb__,mnSz,isTight);
                Bt = [Bt; Btc];
                Lt(:,:,k) = Ltc.*Rc;
                k = k + 1;
            end
            waitbar(i/ns,hw,'searching domains ...');
        end
        close(hw);

        %% find boundaries

        save('structMapDomainsHR','Rdcrop','strInfo','BtbArr','BtbIx','pos'); % high resolution
        save('structMapDomains','Bt','Lt','-v7.3');

        return
        %% display
        RdDisp = RdMb;
        if 2 % display the whole structure
            figure(13);
            %imagesc(Rd); axis image;
            cmn = 64;
            cmnr = 1:cmn;
            CM = jet(cmn);
            [~,ixs] = sort(rand(1,cmn));
            cix = cmnr(ixs);
            CM=CM(cix,:);
            CM(1,:)=0;
            plotBoundaries(BtbArr,RdDisp); 
            colormap(CM); 
            axis image
            set(gcf,'units','pixels','Position',[0,0,size(Rc,2)*mag,size(Rc,1)*mag]); 
            set(gca,'units','pixels','Position',[0,0,size(Rc,2)*mag,size(Rc,1)*mag]);
            imgFig = getframe(gcf);
            dataImg = imgFig.cdata; 
            imwrite(uint16(dataImg),'stDom.tif'); % structMap
        end
        %% struct stats
    end
else % coarse boundaries
    %% pre-process
    %figure(1); clf; % edge selection overlay with rec.

end
%% bleaching times (finds the times the structs bleach)
if 0
    if isempty(fnameRecTimeBlchDIR)
        % select reference region
        refBckGrndFN = 'refBckGrnd.mat';
        if ~exist(refBckGrndFN)
            figure(99)
            nc = 256;
            CM = jet(nc);
            colormap(CM);
            subplot(2,2,1);
            image(A*nc); axis image
            subplot(2,2,2);
            image(Adev*nc); axis image
            title('deviation')
            subplot(2,2,3);
            image(Afamp*nc); axis image
            title('fluct amp')
            subplot(2,2,1);
            maximize;
            hr = imrect(gca);
            p = round(getPosition(hr));
            save(refBckGrndFN,'p')
        else
            load(refBckGrndFN); %p
        end

        % calculate reference intensity values
        nrf = 100; % number of reference frames
        load(['..\fname.mat']); % 
        fname = sprintf('..\\%s',fname); % acq. movie
        fid = fopen('..\frames.txt');
        ftxt = fgetl(fid);
        fr = sscanf(ftxt,'%i:%i');
        ftxt = fgetl(fid);
        acqTime = sscanf(ftxt,'acqTime:%i')/1000;
        nFrBin = round(5/acqTime);
        fr1 = fr(1);
        ns = size(Lt,3);
        hw = waitbar(0,'finding reference intensities');
        fr2 = f2*nFrBin; % f2 from cropTimeWin
        cc = 1;
        for i = fr2-nrf+1:fr2 % find bckgrnd and struct intensity at last frames
            I = double(repelem(imread(fname,i+fr1-1),4,4)); % READ
            Ib = double(I(p(2):p(2)+p(4),p(1):p(1)+p(3))); % bckgrnd
            ib2(cc) = mean(Ib(:));
            for j = 1:ns % each structure
                Is = I.*Lt(:,:,j);
                is(cc,j) = mean(Is(Is>0));
            end
            waitbar(cc/nrf,hw,sprintf('finding reference intensities'));
            cc=cc+1; % counter
        end
        close(hw);
        is = sort(is);
        is = mean(is(1:nrf/2,:),1);
        ib2 = sort(ib2);
        ib2 = mean(ib2(1:nrf/2));
        save('tempBleachRef','is','ib2')

        stBlch = zeros(ns,1);
        wl = 3; % window length [frames]
        Iw = zeros([size(Amean)*4 wl]);
        ibw = zeros(wl,1);
        Frames =10000;
        hw = waitbar(0,'finding bleaching times');
        lastFrame = Frames-fr1+1;
        for i=1:lastFrame % read acq. movie, find bleach times
            I = double(repelem(imread(fname,i+fr1-1),4,4)); % READ
            Ib = double(I(p(2):p(2)+p(4),p(1):p(1)+p(3))); % bckgrnd
            ib = mean(Ib(:));
            ibw = circshift(ibw,1); 
            ibw(1)= ib;
            Iw = circshift(Iw,1,3);
            Iw(:,:,1) = I;
            if i<wl, continue; end
            for j = 1:ns
                if stBlch(j), continue; end; % already bleached
                st = repmat(Lt(:,:,j),[1 1 wl]).*Iw;
                stmx = max(st,[],3);
                thresh(j) = mean(ibw)/ib2*is(j);
                if isempty(find(stmx>thresh(j)))
                    stBlch(j)=1;
                    stBlchFrm(j) = i;
                    waitbar(sum(stBlch)/numel(stBlch),hw,sprintf('finding bleaching times. frame #%i, bleach ratio:%i/%i',i,sum(stBlch),numel(stBlch)))
                end
            end
            if sum(1-stBlch) == 0, break; end; % all bleached at this frame
        end
        close(hw);
        if ~sum(1-stBlch) == 0, find((stBlch==0)); error('doesnt bleach'); end; 
        stNonBlch = find((stBlch==0));
        stBlchFrm(stNonBlch) = lastFrame; 


        % generate and print images
        fnXYZtraceData0 = rdir('..\_**\*traceData0-coeff*');
        load(fnXYZtraceData0.name); % frmNoTrace trInf
        frTr = trInf(:,1);
        stNoImg_ = 1:ns;
        stNoImg = repmat(stNoImg_',[1 size(Lt,2) size(Lt,1)]);
        stNoImg = permute(stNoImg,[2 3 1]);
        stNoImg = stNoImg.*Lt;

        nTr = size(trInf,1); % number of traces
        trImg = zeros(size(Lt(:,:,1))); % trace image (single point)
        stTr = cell(ns,1);
        for i = 1:nTr % all traces
            trImg(trInf(i,10),trInf(i,9)) = 1;
            stNo = find(trImg.*stNoImg>0);
            stTr{stNo} = [stTr{stNo} i]; % traces in the structure
        end

        recSumBlch = zeros(size(trImg));
        recTimeBlch = zeros(size(trImg));
        pxBin = 0.25; % pixel size multiplier
        for i = 1: ns % each structure
            fr1 = stBlchFrm(i);
            fr2 = f2; % from cropTimeWin
            trSelTime = find((fr1<=frTr) .* (frTr<=fr2)); % selected traces (5sec bin)
            trSelXY = stTr{stNo};
            trSel_ = sort([trSelTime trSelXY]);
            trSel = find(histcounts(trSel_,(0:nTr)+0.5)==2);
            recXbin = trInf(trSel,9);
            recYbin = trInf(trSel,10);
            recInt = trInf(trSel,11);
            rX = round(recXbin/pxBin);
            rY = round(recYbin/pxBin);
            rXY = sub2ind(size(trImg),rX,rY);

            [stInt,~,bin] = histcounts(rXY,(0:numel(trImg))+0.5); % intensity (# recs)
            rXY = unique(rXY);
            recSumBlch(rXY) = stInt;
            recTimeBlch(rXY) = accumarray(bin', frTr, [nTr 1], @max, 0);

        end
        recSumBlchFN = ('..\binImgRcrtSum_stBlchTime.tif');
        recTimeBlchFN = ('..\binImgRcrtTime_stBlchTime.tif');
        imwrite(uint16(recSumBlch),recSumBlchFN);
        imwrite(uint16(recTimeBlch),recTimeBlchFN);
    else
        fnameRecSumBlch = fnameRecSumBlchDIR(1).name;
        fnameRecTimeBlch = fnameRecTimeBlchDIR(1).name;
        recSumBlch = imread(fnameRecSumBlch,1);
        recTimeBlch = imread(fnameRecTimeBlch,1);
    end
    %% find coarse boundaries
    R = double(recSumBlch);
    [Bt,Lt] = findBoundaries(R,mnSz0,0);
    save(stMapMATFN,'Lt','Bt','-v7.3') % structMap

    dispBoundaries(R,Bt,Lt,mag,[],isSelSt);
    imgFig = getframe(gcf);
    dataImg = imgFig.cdata; 
    imwrite(uint16(dataImg),stMapImgFN); % structMap

end

%% assign structure types    =====================================
if isAutoFindDomains % assign domains to structure type

    posCtrl2 = [-20,120,120,180];
    while 1 % escape to quit
        figure(5);figure(11);
        % select a structure
        set(gcf,'PointerShapeCData',ones(16)+1)
        [px,py]=ginput(1);

        kp=get(gcf,'CurrentCharacter')+1;
        if kp == 28, break; end; % escape for quit

        px=round(px);
        py=round(py);
        inImg = zeros(size(R));
        inImg(py,px) = 1;
        ixS = []; 
        cS = 0; 
        while isempty(ixS) % continue till hit
            for i=1:ns % each structure
                Lt = double(RdcMb==i);
                if sum(sum(inImg.* conv2(Lt,ones(cS),'same')))
                    ixS = i; 
                    break;
                end
            end
            cS = cS+3;
        end
        Rdc = Rdcrop{ixS};
        p = pos(ixS,:);
        BtbCrop = BtbArr(BtbIx{ixS});
        for j = 1:numel(BtbCrop)
            BtbC = BtbCrop{j};
            BtbC = [BtbC(:,1)-p(1)+1 BtbC(:,2)-p(3)+1];
            BtbCrop{j} = BtbC;
        end


        CM65 = colormap('parula');
        CM65(65,:)= 0;
        colormap(CM65);
 
        colormap('parula')
        Rdc0 = Rdc;
        Rdc = ceil(Rdc/max(Rdc(:))*64);
        figure(1);
        mag=8; %pause(0.1)
        plotBoundaries(BtbCrop,Rdc,0); 
        set(gcf,'units','pixels','Position',[120,120,size(Rdc,2)*mag,size(Rdc,1)*mag]); 
        set(gca,'units','pixels','Position',[0,0,size(Rdc,2)*mag,size(Rdc,1)*mag]);

        CM = colormap;
        cmix = size(CM,1); % # of colormap indices
        Rs = R(p(1):p(2),p(3):p(4)); % cropped section
        As = A(p(1):p(2),p(3):p(4))*cmix; % cropped image
        Afamps = Afamp(p(1):p(2),p(3):p(4))*cmix;
        Adevs = Adev(p(1):p(2),p(3):p(4))*cmix;
        Rtimes = Rtime(p(1):p(2),p(3):p(4));

        figure(4); % Amean overlay
        set(gcf,'name','intensity mean')
        image(As)
        hold on;
        [Y,X] = find(Rs>0);
        scatter(X,Y,'.','r');
        hold off;
        set(gcf,'units','pixels','Position',[120+size(As,2)*mag+10,120,size(As,2)*mag,size(As,1)*mag]); 
        set(gca,'units','pixels','Position',[0,0,size(As,2)*mag,size(As,1)*mag]);

        figure(6); % Afamp overlay
        set(gcf,'name','flicker amplitude')
        image(Afamps)
        hold on;
        scatter(X,Y,'.','r');
        hold off;
        set(gcf,'units','pixels','Position',[120+2*(size(As,2)*mag+10),120,size(As,2)*mag,size(As,1)*mag]); 
        set(gca,'units','pixels','Position',[0,0,size(As,2)*mag,size(As,1)*mag]);   

        figure(8); % Adev overlay
        set(gcf,'name','intensity deviation')
        image(Adevs)
        hold on;
        scatter(X,Y,'.','r');
        hold off;
        set(gcf,'units','pixels','Position',[120+3*(size(As,2)*mag+10),120,size(As,2)*mag,size(As,1)*mag]); 
        set(gca,'units','pixels','Position',[0,0,size(As,2)*mag,size(As,1)*mag]); 

        Rtimes = Rtime(p(1):p(2),p(3):p(4));
        figure(10); % time progress
        set(gcf,'name','intensity deviation')
        image(Rtimes)
        hold on;
        scatter(X,Y,'.','r');
        hold off;
        set(gcf,'units','pixels','Position',[120+4*(size(As,2)*mag+10),120,size(As,2)*mag,size(As,1)*mag]); 
        set(gca,'units','pixels','Position',[0,0,size(As,2)*mag,size(As,1)*mag]); 

        figure(3);
        imagesc(Rs);
        set(gcf,'units','pixels','Position',[120,120,size(Rdc,2)*mag,size(Rdc,1)*mag]); 
        set(gca,'units','pixels','Position',[0,0,size(Rdc,2)*mag,size(Rdc,1)*mag]);

        for i=1:ns
            RdcDisp = Rdc;
            RdcDisp(Rdc0==i)=65;
            figure(1); clf; % edge selection overlay with rec.
            plotBoundaries(BtbCrop,RdcDisp,0); 
            set(gcf,'units','pixels','Position',[120,120,size(Rdc,2)*mag,size(Rdc,1)*mag]); 
            set(gca,'units','pixels','Position',[0,0,size(Rdc,2)*mag,size(Rdc,1)*mag]);

            colormap(CM65);
            delete(5);
            figure(5) % select structure type                      
            ss = {'1: str. body',...
                '2: str. growth',...
                '3: str. pit',...
                '4: str. hotspot',...
                '5: pit',...
                '6: hotspot'};
            set(gcf,'units','pixels','Position',posCtrl2);
            t6 = uicontrol('Style','text',...
                       'Position',[30 120 210 20],...
                       'HorizontalAlignment','left',...
                       'String',ss{6});
            t5 = uicontrol('Style','text',...
                       'Position',[30 100 210 20],...
                       'HorizontalAlignment','left',...
                       'String',ss{5});
            t4 = uicontrol('Style','text',...
                       'Position',[30 80 210 20],...
                       'HorizontalAlignment','left',...
                       'String',ss{4});
            t3 = uicontrol('Style','text',...
                       'Position',[30 60 210 20],...
                       'HorizontalAlignment','left',...
                       'String',ss{3});
            t2 = uicontrol('Style','text',...
                       'Position',[30 40 210 20],...
                       'HorizontalAlignment','left',...
                       'String',ss{2});
            t1 = uicontrol('Style','text',...
                       'Position',[30 20 210 20],...
                       'HorizontalAlignment','left',...
                       'String',ss{1});
            % wait for selection (key)       
            ks = []; pause(0.1);
            figure(4); figure(2);figure(5)
            while isempty(ks), pause(0.1); ks=get(gcf,'CurrentCharacter'); end; 
            figure(5)
            if ks=='s' % save figures
                figs = [1 3 4 6 8];
                imgFN=sprintf('stImg%03i_%iX%iY%ix%i.tif',ixS,p(1),p(3),p(2)-p(1)+1,p(4)-p(3)+1);
                if exist(imgFN), delete(imgFN); end
                figure(1);
                plotBoundaries(BtbCrop,Rdc,0); 
                for i = 1:numel(figs)
                    figure(figs(i));
                    imgFig = getframe(gcf);
                    imgOut = imgFig.cdata;
                    imwrite(imgOut,imgFN,'WriteMode','append','Compression', 'none') 
                end
                break
            else
                set(eval(sprintf('t%i',str2num(ks))),'BackgroundColor',[1 1 0]); % highlight text
                eval(sprintf('ps(end-ns+%i,2)=%i',i,str2num(ks))); % register structure type
                %ps_(:,2:end) = ps(end-ns+1:end,2:end);
            end


        end            

%(Rdcrop strInfo BtbArr BtbIx pos)

    end

else % coarse boundaries

    %% =========================================        
    %% RE-MAP (R,Bt,Lt) (interactive mode)
    % select stuctures 
    LtNew = Lt; 
    BtNew = Bt; clear Bt Lt;
    if isempty(ps), 
        ps = [1:size(LtNew,3)]'; 
        ps(:,3)=ps(:,1); 
        %ps(:,4)=f1; 
        %ps(:,5)=f2; 
    end
    Rnew = R;
    mag=8;
    if isSelSt % indices of processed structures
        [ixReMap,~]=find((ps(:,2)==1) + (ps(:,2)==11));
        load(edgeDistDataDIR(1).name);
        rsel=0;
        fid = fopen('edgeDistMeanArr.txt','wt');
        edgeDistMeanArr_ = edgeDistMeanArr;
        edgeDistMeanArr_(:,end+1) = 1:size(edgeDistMeanArr,1);
        edgeDistMeanArr_(:,end+1) = edgeDistMeanArr_(:,1)./edgeDistMeanArr_(:,2);
        edgeDistMeanArr_ = circshift(edgeDistMeanArr_,2,2);
        fprintf(fid,'%i %.02f %.02f %.02f %.02f %i\n',edgeDistMeanArr_');
        fclose(fid);
    end
    while 1 % escape to quit
        Ns=size(LtNew,3); % number of structures
        figure(888)
        figure(1)

        % select a structure
        set(gcf,'PointerShapeCData',ones(16)+1);
        isBreak = 0;
        while 1
            [px,py,button]=ginput(1);

            kp=get(gcf,'CurrentCharacter')+1;
            if kp == 28, isBreak = 1; break; end; % escape for quit
            if button
                break
            end
        end
        if isBreak, break; end;

        px=round(px);
        py=round(py);
        inImg = zeros(size(R));
        inImg(py,px) = 1;
        ixS = []; 
        cS = 0; 
        while isempty(ixS) % continue till hit
            for i=1:Ns % each structure
                if sum(sum(inImg.* conv2(LtNew(:,:,i),ones(cS),'same')))
                    ixS = i; 
                    break;
                end
            end
            cS = cS+3;
        end


        % display structure
        b = BtNew{ixS};
        fs = 10; % frame size
        y1 = min(b(:,1))-fs; if y1<1, y1=1; end
        y2 = max(b(:,1))+fs; if y2>size(R,1), y2=size(R,1); end
        x1 = min(b(:,2))-fs; if x1<1, x1=1; end
        x2 = max(b(:,2))+fs; if x2>size(R,2), x2=size(R,2); end
        Ltc = LtNew(:,:,ixS);
        Rc = R.*Ltc;
        Rc = Rc(y1:y2,x1:x2); % cropped struct 
        Rs = R(y1:y2,x1:x2); % cropped section
        CM = colormap;
        cmix = size(CM,1); % # of colormap indices
        As = A(y1:y2,x1:x2)*cmix; % cropped image
        Afamps = Afamp(y1:y2,x1:x2)*cmix;
        Adevs = Adev(y1:y2,x1:x2)*cmix;
        Rtimes = Rtime(y1:y2,x1:x2);

        
    
    
        figure(4); % Amean overlay
        set(gcf,'name','intensity mean')
        image(As); axis image;
        hold on;
        [Y,X] = find(Rs>0);
        scatter(X,Y,'.','r');
        hold off;
        set(gcf,'units','pixels','Position',[120+size(As,2)*mag+10,120,size(As,2)*mag,size(As,1)*mag]); 
        set(gca,'units','pixels','Position',[0,0,size(As,2)*mag,size(As,1)*mag]);             

        figure(6); % Afamp overlay
        set(gcf,'name','flicker amplitude')
        image(Afamps); axis image;
        hold on;
        scatter(X,Y,'.','r');
        hold off;
        set(gcf,'units','pixels','Position',[120+2*(size(As,2)*mag+10),120,size(As,2)*mag,size(As,1)*mag]); 
        set(gca,'units','pixels','Position',[0,0,size(As,2)*mag,size(As,1)*mag]);   

        figure(8); % Adev overlay
        set(gcf,'name','intensity deviation')
        image(Adevs); axis image;
        hold on;
        scatter(X,Y,'.','r');
        hold off;
        set(gcf,'units','pixels','Position',[120+3*(size(As,2)*mag+10),120,size(As,2)*mag,size(As,1)*mag]); 
        set(gca,'units','pixels','Position',[0,0,size(As,2)*mag,size(As,1)*mag]); 

        figure(10); % time progress
        set(gcf,'name','recruitment times')
        image(Rtimes)
        set(gcf,'units','pixels','Position',[120+4*(size(As,2)*mag+10),120,size(As,2)*mag,size(As,1)*mag]); 
        set(gca,'units','pixels','Position',[0,0,size(As,2)*mag,size(As,1)*mag]); 
        colormap('jet')
        

        rsel = 1; % default method#1 (remove recs) % press ESC to quit
        if isSelSt % indices of processed structures
            %[ixReMap,~]=find((ps(:,2)==1) + (ps(:,2)==11));
            [ixReMap,~]=find((ps(:,2)==11));
            load(edgeDistDataDIR(1).name);
            rsel=0;
        end

        figure(3);figure(5);

        posCtrl = [-20,120,250,240];
        posCtrl2 = [-20,120,120,180];
        nsr = 0; % number of sub-regions
        pbsel = 0;
        isTightenBound = 0;
        tightPosNew = tightPos{ixS};
        tightPosAdd = [];
        % re-map the selected structure
        while 1 % press esc or enter to quit
            Btn = {};
            Ltn = [];
            i=0;
            tightPosNew = [tightPosNew ; tightPosAdd];
            while i < size(Rc,3) % each sub-structure
                if isAutoBorderDetection
                    b_ = b;
                    b_(:,1) = b_(:,1)-y1+1;
                    b_(:,2) = b_(:,2)-x1+1;
                    Btn = {b_};
                    Ltn = Ltc(y1:y2,x1:x2);
                    break
                end
                i=i+1;
                [Btn_,Ltn_] = findBoundaries(Rc(:,:,i),mnSz,0); % new structures
                tps = tightPosNew;
                if ~isempty(tps)
                    tps(:,1) = tps(:,1)-x1+1;
                    tps(:,2) = tps(:,2)-y1+1;
                    [Btn_,Ltn_] = tightenBoundaries(Btn_,Ltn_,Rs,tps);
                end
                if isempty(Ltn_), Rc(:,:,i) = []; continue; end;
                ns2 = numel(Btn_);
                %for i = 1:ns2
                    Btn = [Btn;Btn_];
                    if isempty(Ltn)
                        Ltn = Ltn_;
                    else
                        Ltn(:,:,end+1:end+ns2) = Ltn_;
                    end
                %end
            end
            
            figure(999)
            set(gcf,'units','pixels','Position',[320+3*(size(As,2)*mag+10),size(As,1)*mag+150,size(As,2)*mag,size(As,1)*mag]); 
            set(gca,'units','pixels','Position',[0,0,size(As,2)*mag,size(As,1)*mag]);         
            %imagesc(sum(Ltn,3))
            imagesc(sum(Ltc,3))

            figure(2); clf; % edge selection overlay with rec.
            dispBoundaries(Rs,Btn,Ltn,mag,[],isSelSt)
            if isSelSt % display specs
                delete(3)
                figure(3) % controls
                set(gcf,'units','pixels','Position',posCtrl);
                
                [ixEdgeProc,~]=find(ixReMap==ixS);
                btn1 = uicontrol('Style','pushbutton','String','remove structure','Position',[20 5 130 20],...
                          'HandleVisibility','off','Callback','bsel=11;','HorizontalAlignment','left');
                yy = edgeDistRec(ixEdgeProc,:);
                xx = 0:size(edgeDistRef,2)-1;
                szdf=size(edgeDistRef,2)-size(edgeDistRec,2);
                if szdf, yy(end+1:end+szdf)=0; end;
                plot(xx,yy,'b');
                hold on; 
                plot(xx,edgeDistRef(ixEdgeProc,:),'r'); 
                plot(xx,edgeDistStr(ixEdgeProc,:),'k'); 
                hold off;
                title(sprintf('rec=%.02f int=%.02f; uni=%.02f',edgeDistMeanArr(ixEdgeProc,[1 3 2])*160/1.5/4))
                set(gca,'units','pixels','Position',posCtrl+[65 -75 -60 -70]);   
                grid minor
               
           else % re-map structure
                delete(3)
                figure(3) % controls
                set(gcf,'units','pixels','Position',posCtrl);
                dt = 20;
                btn1 = uicontrol('Style','pushbutton','String','remove recs (#1)','Position',[20 150+dt 130 20],...
                          'HandleVisibility','off','Callback','bsel=1;','HorizontalAlignment','left');
                btn2 = uicontrol('Style','pushbutton','String','select sub-struct (#2)','Position',[20 120+dt 130 20],...
                          'HandleVisibility','off','Callback','bsel=2;','HorizontalAlignment','left');
                btn3 = uicontrol('Style','pushbutton','String','tighten the boundary','Position',[20 90+dt 130 20],...
                          'HandleVisibility','off','Callback','bsel=5;','HorizontalAlignment','left');
                btn4 = uicontrol('Style','pushbutton','String','finish','Position',[20 60+dt 130 20],...
                          'HandleVisibility','off','Callback','bsel=3;','HorizontalAlignment','left');
                btn5 = uicontrol('Style','pushbutton','String','cancel','Position',[20 30+dt 130 20],...
                          'HandleVisibility','off','Callback','bsel=4;','HorizontalAlignment','left');
                txt = uicontrol('Style','text','String','finish to save','Position',[20 -10+dt 130 40]);
            end
            
            if ~isempty(tightPos{ixS})
                if exist('btn2'), delete(btn2); end;
            end
            bsel = 0;
            kp=[]; % press anykey 
            %set(gcf,'CurrentCharacter',char(1))
            pause(0.5)
            if ~rsel, while ~bsel && isempty(kp), pause(0.1); kp=get(gcf,'CurrentCharacter')+1; end; elseif ~isSelSt, bsel = 1; end;
            if kp == 28, bsel=4; elseif kp == 14, bsel =3; end; % esc and enter keys
            if (kp == 14) & logical(isSelSt), bsel = 11; end; % enter key in isSelSt mode
            if bsel == 11 % isSelSt:1 -> remove selected structure
                if ps(ixS,2) == 0, continue; end;
                ps(ixS,2) = 11;
                break;
            elseif bsel==1 %method#1
                rsel = 1; % default method#1 (remove recs) % press ESC to quit
                set(txt,'String','#1 active. Press escape to return');
                figure(2);
                h=imrect(gca);
                if isempty(h)
                    rsel = 0;
                else
                    p = round(getPosition(h));
                    p4=p(2)+p(4);
                    p3=p(1)+p(3);
                    if p(2)+p(4)> size(R,1), p4=size(R,1); end
                    if p(1)+p(3)> size(R,2), p3=size(R,2); end
                    if p(1)<=0 || p(2)<=0, p(1)=0;p(2)=0; end  % BREAK
                    Rc(p(2):p4,p(1):p3,:) = 0;
                    Rs(p(2):p4,p(1):p3,:) = 0;
                    Rnew([p(2):p4]+y1-1,[p(1):p3]+x1-1) = 0;
                end
            elseif bsel == 2 % method#2
                if pbsel == 1
                    set(txt,'String','#2 done before. save or remove recs');
                else
                    pbsel = 1;
                    set(txt,'String','#2 active. Press escape to return');
                    figure(2);
                    BW = roipoly;
                    if isempty(BW), continue; end;
                    RcOld=Rc;
                    Rc(:,:,2) = RcOld.*BW;
                    Rc(:,:,1) = RcOld.*-(BW-1);
                    clear RcOld;
                    rsel = 1; % default method#1 (remove recs) % press ESC to quit
                end
            elseif bsel == 5 % method#3
                ccc = 3;
                isTightenBound = 1;
                while (1)
                    delete(2)
                    figure(2); clf; % edge selection overlay with rec.
                    dispBoundaries(Rs,Btn,Ltn,mag,[],isSelSt)
                    [px,py]=ginput(1);
                    px = round(px);py= round(py);
                    kp=get(gcf,'CurrentCharacter')+1;
                    if kp == 28, break; end; % escape for quit
                    tightPosAdd = [tightPosAdd; px+x1-1 py+y1-1];
                    
                    [Btn,Ltn] = tightenBoundaries(Btn,Ltn,Rs,[px py]);
                    figure(2)
                    dispBoundaries(Rs,Btn,Ltn,mag,[],isSelSt)
                end
                cccc=2;
            elseif bsel == 4 %: cancel
                Btn = {};
                Ltn = [];
                break
            elseif bsel == 3 %: finish struct
                % save
                tightPos{ixS} = tightPosNew;
                ns = size(Ltn,3);
                % switch to global coordinates
                LtN = zeros(size(R,1),size(R,2),ns);
                BtN = [];
                for i = 1:ns % each sub struct
                    BtN{i} = Btn{i}+repmat([y1-1 x1-1],size(Btn{i},1),1);
                    LtN(y1:y2,x1:x2,i) = Ltn(:,:,i);
                end 
                BtNew = BtNew(~ismember(1:size(BtNew,1), ixS)); % structures passed without change 
                BtNew = [BtNew; BtN'];
                tightPos_ = tightPos(ixS);
                tightPos(ixS:end-1)=tightPos(ixS+1:end);
                tightPos(end:end+ns-1)=tightPos_;
                LtNew(:,:,ixS) = [];
                LtNew(:,:,end+1:end+ns)= LtN;
                %LtNew = cat(3,LtNew(:,:,[1:ixS-1 ixS+1:end]),LtN); %memory inefficient
                ns2 = size(LtNew,3);

                ixAdd = ns2-ns+1:ns2;
                nps = size(ps,1);
                %ps(nps:end+nsAdd,1:3) = nan;
                [psRix,~] = find(ixS==ps(1:nps,1));
                if size(ps,1)>psRix, ps(psRix+1:end,1) = ps(psRix+1:end,1)-1; end;
                ixg = ps(psRix,3); % group index;
                ps(psRix,:) = []; % remove
                ps(end+1:end+ns,1) = ixAdd;
                ps(end-ns+1:end,3) = ixg;

                if ps(end,1)-ps(1,1)+1 ~= numel(ps)
                    insertBreakPoint
                    %input('error:ps')
                end

                % register type
                for i=1:ns

                    ps_ = zeros(ns,size(ps,2));
                    ps_(:,1) = [1:ns]; 
                    ps_(:,2:end) = ps(end-ns+1:end,2:size(ps,2));
                    ps_(i,2) = ps_(i,2)+20; % make current struct bold

                    figure(2); clf; % edge selection overlay with rec.
                    dispBoundaries(Rs,Btn,Ltn,mag,ps_,isSelSt);



                    delete(5);
                    figure(5) % select structure type                      
                    ss = {'1: str. body',...
                        '2: str. growth',...
                        '3: growing str.        ',...
                        '4: str. pit',...
                        '5: str. hotspot',...
                        '6: pit',...
                        '7: hotspot'};
                    set(gcf,'units','pixels','Position',posCtrl2);
                    t7 = uicontrol('Style','text',...
                               'Position',[30 120 210 20],...
                               'HorizontalAlignment','left',...
                               'String',ss{7});
                    t6 = uicontrol('Style','text',...
                               'Position',[30 120 210 20],...
                               'HorizontalAlignment','left',...
                               'String',ss{6});
                    t5 = uicontrol('Style','text',...
                               'Position',[30 100 210 20],...
                               'HorizontalAlignment','left',...
                               'String',ss{5});
                    t4 = uicontrol('Style','text',...
                               'Position',[30 80 210 20],...
                               'HorizontalAlignment','left',...
                               'String',ss{4});
                    t3 = uicontrol('Style','text',...
                               'Position',[30 60 210 20],...
                               'HorizontalAlignment','left',...
                               'String',ss{3});
                    t2 = uicontrol('Style','text',...
                               'Position',[30 40 210 20],...
                               'HorizontalAlignment','left',...
                               'String',ss{2});
                    t1 = uicontrol('Style','text',...
                               'Position',[30 20 210 20],...
                               'HorizontalAlignment','left',...
                               'String',ss{1});
                    % wait for selection (key)       
                    ks = []; pause(0.15);
                    figure(4); figure(2);figure(10);figure(5);
                    ksvals = [50 51 52 53 54 55 28]; % [1 2 3 4 5 6 esc]
                    while isempty(ks) || ~ismember(ks+1,ksvals), pause(0.1); ks=get(gcf,'CurrentCharacter'); end; 
                    figure(5)
                    if ks+1 == 28, 
                        ps_(:,2:end) = ps(end-ns+1:end,2:end);
                        break; 
                    end
                    set(eval(sprintf('t%i',str2num(ks))),'BackgroundColor',[1 1 0]); % highlight text
                    eval(sprintf('ps(end-ns+%i,2)=%i;',i,str2num(ks))); % register structure type
                    ps_(:,2:end) = ps(end-ns+1:end,2:end);

                end
                figure(2); clf; % edge selection overlay with rec.
                dispBoundaries(Rs,Btn,Ltn,mag,ps_,isSelSt);
                pause


                break
            end     % process method
        end        % re-map the selected structure



        % display updated mapping
        figure(1); clf; % edge selection overlay with rec.
        dispBoundaries(R,BtNew,LtNew,1,ps,isSelSt);
    end % select structure
ps(:,2)=11; % all structures are OK
    save(stMapUpdMATFN,'BtNew','LtNew','ps','tightPos','-v7.3'); % structMap

    if strcmp(input('quit(enter) or save img(y)?','s'),'y')
        dispBoundaries(R,BtNew,LtNew,1,ps,isSelSt);
        figure(1)
        imgFig = getframe(gcf);
        dataImg = imgFig.cdata; 
        if isSelSt
            imwrite(uint16(dataImg),stMapUpdImgFN2); % structMap (selected structures)
        else
            imwrite(uint16(dataImg),stMapUpdImgFN); % structMap
        end
    end

    disp(cellLabel);

    %pause
end
    
return;  
    
%% print structures
delete(btn1)
FN = sprintf('%s_str%i',cellLabel,ixS);
figure(2)
        imgFig = getframe(gcf);
        dataImg = imgFig.cdata; 
        imwrite(uint16(dataImg),[FN '_boundary.tif']);
figure(4)
        imgFig = getframe(gcf);
        dataImg = imgFig.cdata; 
        imwrite(uint16(dataImg),[FN '_overlap.tif']);
figure(3)
        imgFig = getframe(gcf);
        dataImg = imgFig.cdata; 
        imwrite(uint16(dataImg),[FN '_hist.tif']);
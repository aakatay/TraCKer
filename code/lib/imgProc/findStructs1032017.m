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
selMask = 0;
%selMask = 0; % no cropping
mnSz0 = 65; % only large structures are processed
mnSz = mnSz0; % only large structures are processed
%mnSz0 = 20; % only large structures are processed
mnDensity = 0.25; % only dense recruiting structs are processed
mnDensity = 0.1;

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
stMapUpdEdgeImgFN= [ 'structMapUpd-Edgeness' cellLabel '.tif']; % map
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
    mag = 1;
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

if selMask | exist('posMask.mat')
    load posMask
    A = A.*double(mask);
    R = R.*double(mask);
end


%% =========== find boundaries =============
mag=1;
if isCropTimeWin | (~exist(stMapMATFN)  & ~exist(stMapUpdMATFN)) % find coarse boundaries
    %R = double(im2bw(R,0));
    mnFilt = [mnSz0 mnDensity];
    [Bt,Lt] = findBoundaries(R,mnFilt);
    
    tightPos = cell(size(Lt,3),1); % selected tight requests
    save(stMapMATFN,'Lt','Bt','tightPos','-v7.3') % structMap

    dispBoundaries(R,Bt,Lt,mag,[],isSelSt);
    figure(1);
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

else % coarse boundaries
    %% pre-process
    %figure(1); clf; % edge selection overlay with rec.

end
%% bleaching times (finds the times the structs bleach)


%% assign structure types    =====================================
if isAutoFindDomains % assign domains to structure type

else % coarse boundaries

    %% =========================================        
    %% RE-MAP (R,Bt,Lt) (interactive mode)
    % select stuctures 
    LtNew = Lt; 
    BtNew = Bt; clear Bt Lt;
    if isempty(ps), 
        ps = [1:size(LtNew,3)]'; 
        ps(:,3) = ps(:,1); 
        ps(:,4) = f1; % time window
        ps(:,5) = f2; 
    elseif size(ps,2)<5
        ps(:,4) = f1; % time window
        ps(:,5) = f2; 
    end
    R = cropTimeWinStructure(binFN,ps,LtNew);
    %Rnew = R;
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

        posCtrl = [-20,120,250,280];
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
                if 0&isAutoBorderDetection
                    b_ = b;
                    b_(:,1) = b_(:,1)-y1+1;
                    b_(:,2) = b_(:,2)-x1+1;
                    Btn = {b_};
                    Ltn = Ltc(y1:y2,x1:x2);
                    break
                end
                i=i+1;
                [Btn_,Ltn_] = findBoundaries(Rc(:,:,i),mnSz); % new structures
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
            
            imagesc(sum(Ltn,3))

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
                btn6 = uicontrol('Style','pushbutton','String','timeSel','Position',[20 180+dt 130 20],...
                          'HandleVisibility','off','Callback','bsel=9;','HorizontalAlignment','left');
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
            if bsel == 9 % time crop
                cP = [x1 y1 x2 y2]; % crop pos
                [f1,f2,binImgLastTime,binImgCr] = cropTimeWinCore(binFN,f1,f2,cP);
                ps(ixS,4:5) = [f1 f2]; % move the rest to (bsel == 3 %: finish struct)
                R(cP(2):cP(4),cP(1):cP(3)) = binImgCr; 
                Rc = binImgCr;
                Rs = binImgCr;
            elseif bsel == 11 % isSelSt:1 -> remove selected structure
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
                    %Rnew([p(2):p4]+y1-1,[p(1):p3]+x1-1) = 0;
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
                f1_ = ps(psRix,4); 
                f2_ = ps(psRix,5); 
                ps(psRix,:) = []; % remove
                ps(end+1:end+ns,1) = ixAdd;
                ps(end-ns+1:end,3) = ixg;
                ps(end-ns+1:end,4) = f1_;
                ps(end-ns+1:end,5) = f2_;

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

    %if strcmp(input('quit(enter) or save img(y)?','s'),'y')
        dispBoundaries(R,BtNew,LtNew,1,ps,isSelSt);
        figure(1)
        imgFig = getframe(gcf);
        dataImg = imgFig.cdata; 
        if isSelSt
            imwrite(uint16(dataImg),stMapUpdImgFN2); % structMap (selected structures)
        else
            imwrite(uint16(dataImg),stMapUpdImgFN); % structMap
        end
        fh=findobj(0,'Number',888);
        if ~isempty(fh) % edgeness
            figure(888)
            imgFig = getframe(gcf);
            dataImg = imgFig.cdata; 
            imwrite(uint16(dataImg),stMapUpdEdgeImgFN);
        end
            
        
    %end

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
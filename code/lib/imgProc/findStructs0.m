% run in recComp\ folder
%% 1- select a rectangle to remove outside dots (stay in image size)
%% 2- to quit select a rect towards outside on upper left of the image 
%% 3- it will find edge and distance to the edge with reference and plot and write to xls
clear
%close all
selMask = 0; % no cropping
selMask = -1; % negative mask
selMask = 1;
mnSz0 = 50; % only large structures are processed
mnSz = 0; % only large structures are processed
fnameRecSumDIR = rdir('..\*binImgRcrtSum_time*');
fnameADIR = rdir('..\*AVG_lap*'); % pre bleach image
fnameMask = 'posMask.mat';
A = imread(fnameADIR.name);
A = repelem(A,4,4);

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

%%
ps= [];
for j = 1:numel(fnameRecSumDIR) % each cell
    close all;
    fnameRecSum = fnameRecSumDIR(j).name;
    stMapImgFN= [ 'structMap' cellLabel '.tif']; % map
    stMapUpdImgFN= [ 'structMapUpd' cellLabel '.tif']; % map
    stMapMATFN= [ 'structMap' cellLabel '.mat']; % map
    stMapUpdMATFN= [ 'structMapUpd' cellLabel '.mat']; % map
    overlayFullFN= [ 'overlayFull' cellLabel '.tif']; % map
    R = double(imread(fnameRecSum));
    
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
    Amag = double(repelem(A,mag,mag));
    Amag = round(Amag*Clevel/max(Amag(:)));
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
    imwrite(AcolImg,overlayFullFN,'Compression','none')
    
    %%  define a mask to process in multiple steps
    if ~exist(fnameMask)
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
    %%
    figure(1); clf; % edge selection overlay with rec.
    mag=1;
    if ~exist(stMapMATFN)  && ~exist(stMapUpdMATFN)
        %R = double(im2bw(R,0));
        [Bt,Lt] = findBoundaries(R,mnSz0);
        save(stMapMATFN,'Lt','Bt','-v7.3') % structMap
        
        dispBoundaries(R,Bt,Lt,mag,[]);
        imgFig = getframe(gcf);
        dataImg = imgFig.cdata; 
        imwrite(uint16(dataImg),stMapImgFN); % structMap
    elseif exist(stMapUpdMATFN) % load updated file
        load(stMapUpdMATFN);
        Lt = LtNew;
        Bt = BtNew; clear BtNew LtNew;
        dispBoundaries(R,Bt,Lt,mag,ps);
    else % exist(stMapMATFN)
        load(stMapMATFN);
        dispBoundaries(R,Bt,Lt,mag,ps);
    end % already mapped once

    
    
    %% RE-MAP (R,Bt,Lt)
    % select stuctures
    LtNew = Lt;
    BtNew = Bt; clear Bt Lt;
    Rnew = R;
    mag=8;
    while 1 % click empty region to quit
        Ns=size(LtNew,3); % number of structures
        figure(1)
        
        % select a structure
        set(gcf,'PointerShapeCData',ones(16)+1)
        [px,py]=ginput(1);
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
        As = A(y1:y2,x1:x2); % cropped image
        
        figure(4)
        imagesc(As)
        hold on;
        [Y,X] = find(Rs>0);
        scatter(X,Y,'.','r');
        hold off;
        set(gcf,'units','pixels','Position',[120+size(As,2)*mag+10,120,size(As,2)*mag,size(As,1)*mag]); 
        set(gca,'units','pixels','Position',[0,0,size(As,2)*mag,size(As,1)*mag]);             



        figure(3);figure(5);
        
        posCtrl = [-20,120,250,150];
        rsel = 1; % default method#1 (remove recs) % press ESC to quit
        nsr = 0; % number of sub-regions
        % re-map the selected structure
        while 1 
            Btn = {};
            Ltn = [];
            i=0;
            while i < size(Rc,3) % each sub-structure
                i=i+1;
                [Btn_,Ltn_] = findBoundaries(Rc(:,:,i),mnSz); % new structures
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
            
            figure(2); clf; % edge selection overlay with rec.
            dispBoundaries(Rs,Btn,Ltn,mag,[])

            delete(3)
            figure(3) % controls
            set(gcf,'units','pixels','Position',posCtrl);
            btn1 = uicontrol('Style','pushbutton','String','remove recs','Position',[20 120 160 20],...
                      'HandleVisibility','off','Callback','bsel=1;');
            btn2 = uicontrol('Style','pushbutton','String','select sub-struct','Position',[20 90 160 20],...
                      'HandleVisibility','off','Callback','bsel=2;');
            btn3 = uicontrol('Style','pushbutton','String','finish','Position',[20 60 160 20],...
                      'HandleVisibility','off','Callback','bsel=3;');
            btn3 = uicontrol('Style','pushbutton','String','cancel','Position',[20 30 160 20],...
                      'HandleVisibility','off','Callback','bsel=4;');
                  
            bsel = 0;
            kp=[]; % press anykey 
            %set(gcf,'CurrentCharacter',char(1))
            pause(0.5)
            if ~rsel, while ~bsel && isempty(kp), pause(0.1); kp=get(gcf,'CurrentCharacter')+1; end; else, bsel = 1; end;
            if bsel==1 %method#1
                rsel = 1; % default method#1 (remove recs) % press ESC to quit
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
                figure(2);
                BW = roipoly;
                if isempty(BW), continue; end;
                RcOld=Rc;
                Rc(:,:,2) = RcOld.*BW;
                Rc(:,:,1) = RcOld.*-(BW-1);
                clear RcOld;
                rsel = 1; % default method#1 (remove recs) % press ESC to quit
            elseif bsel == 4 %: cancel
                Btn = {};
                Ltn = [];
                break
            else %bsel == 3 %: finish struct
                % save
                ns = size(Ltn,3);
                % switch to global coordinates
                LtN = zeros(size(R,1),size(R,2),ns);
                for i = 1:ns % each sub struct
                    BtN{i} = Btn{i}+repmat([y1-1 x1-1],size(Btn{i},1),1);
                    LtN(y1:y2,x1:x2,i) = Ltn(:,:,i);
                end 
                BtNew = BtNew(~ismember(1:size(BtNew,1), ixS)); % structures passed without change 
                BtNew = [BtNew; BtN'];
                LtNew = cat(3,LtNew(:,:,[1:ixS-1 ixS+1:end]),LtN);
                ns2 = size(LtNew,3);
                
                ixAdd = ns2-ns+1:ns2;
                nps = size(ps,1);
                %ps(nps:end+nsAdd,1:3) = nan;
                [psRix,~] = find(ixS==ps(1:nps));
                if isempty(psRix) % first time selected
                    if ~isempty(ps), ps(:,1) = ps(:,1)-1; end;
                    ps(end+1:end+ns,1) = ixAdd;
                else % re-selected
                    ps0 = ps(1:psRix-1,:);
                    ps2 = [[ps(psRix+1:nps,1)-1]' ps(psRix+1:nps,2:end)'];
                    ps = [ps0 ps2]; % remove
                    ps(:,1) = [ps(1:end,1) ixAdd]';
                end
                if size(ps,2)<3, ps(:,2:3) = nan; end
                
                if ps(end,1)-ps(1,1)+1 ~= numel(ps)
                    insertBreakPoint
                    %input('error:ps')
                end

                % register type
                for i=1:ns
                    
                    ps_ = nan(ns,size(ps,2));
                    ps_(:,1) = [1:ns]; 
                    ps_(:,2:end) = ps(end-ns+1:end,2:size(ps,2));
                    ps_(i,2) = ps_(i,2)+10; % make current struct bold
                    
                    figure(2); clf; % edge selection overlay with rec.
                    dispBoundaries(Rs,Btn,Ltn,mag,[])

                    delete(5);
                    figure(5) % select structure type
                    set(gcf,'units','pixels','Position',posCtrl);
                    t1 = uicontrol('Style','text',...
                               'Position',[30 100 210 20],...
                               'HorizontalAlignment','left',...
                               'String','1: main body');
                    t2 = uicontrol('Style','text',...
                               'Position',[30 80 210 20],...
                               'HorizontalAlignment','left',...
                               'String','2: grow: body');
                    t3 = uicontrol('Style','text',...
                               'Position',[30 60 210 20],...
                               'HorizontalAlignment','left',...
                               'String','3: grow: edge');
                    t4 = uicontrol('Style','text',...
                               'Position',[30 40 210 20],...
                               'HorizontalAlignment','left',...
                               'String','4: grow: branch');
                    t5 = uicontrol('Style','text',...
                               'Position',[30 20 210 20],...
                               'HorizontalAlignment','left',...
                               'String','5: btw bodies');
                    % wait for selection (key)       
                    ks = []; pause(0.1);
                    while isempty(ks), pause(0.1); ks=get(gcf,'CurrentCharacter'); end; 
                    figure(5)

                    set(eval(sprintf('t%i',str2num(ks))),'BackgroundColor',[1 1 0]); % highlight text
                    eval(sprintf('ps(end-ns+%i,2)=%i',i,str2num(ks))); % register structure type
                    
                    ps_(:,2:end) = ps(end-ns+1:end,2:end);
                end
                figure(2); clf; % edge selection overlay with rec.
                dispBoundaries(Rs,Btn,Ltn,mag,ps_)
                pause
                        
                        
                break
            end
        end        
        

        
        % display updated mapping
        figure(1); clf; % edge selection overlay with rec.
        dispBoundaries(R,BtNew,LtNew,1,ps);
    end
    
    save(stMapUpdMATFN,'BtNew','LtNew','ps') % structMap
    
    if strcmp(input('quit(enter) or save img(y)?','s'),'y')
        dispBoundaries(R,BtNew,LtNew,1,[]);
        imgFig = getframe(gcf);
        dataImg = imgFig.cdata; 
        imwrite(uint16(dataImg),stMapUpdImgFN); % structMap
    end

    disp(cellLabel);
    
    %pause
    
end

    close all;
    clear
    F = findall(0,'type','figure','tag','TMWWaitbar'); delete(F)
    addpath('C:\MATLAB\TraCKer\code\strainCode\Kirchausen')    
    
    inputInfo = dir('inputInfo*.mat');
    load(inputInfo.name);
    dx = PixelSize; 
    dy = dx;
    dz = PlaneDist;

    isCropXY = 0;
    isCombTraces = 1; 
    
        
    sptJmp = 6; %[pixels] % the min. distance between spots of a trace in consecutive frames (spot jump)
    sptReAppearTime = 2; %[frames] 
    minTraceLength = 4; % [frames] traces shorter than this value are not plotted
    
    isPixelsNonIsotropic = 0;
    
    if isPixelsNonIsotropic
        getDeSkewedPixelSizes = 1;
        DeSkewImage; %getDeSkewedPixelSizes (dx, dy, dz)
        WindowSizeX = 3;
        WindowSizeY = 5;
        BigWindowSizeX = 5;
        BigWindowSizeY = 9;
    else
        WindowSizeX = WindowSize;
        WindowSizeY = WindowSizeX;
        BigWindowSizeX = WindowSize + 5;
        BigWindowSizeY = BigWindowSizeX;
    end
    
    

    %READ the ORIGINAL file
    if exist('fname')
        ;
    elseif exist('fname.mat')
        load fname
        if ~exist(fname)
            fname = sprintf('..\\%s',fname);
        end
    else
        fname = 'data.tif';
        if ~exist(fname), fname2 = rdir(sprintf([fname(1:numel(fname)-4) '*.tif'])); fname=fname2.name; end;
        save('fname','fname')
    end
    imgFrst = imread(fname);
    [Boy1,En1]=size(imgFrst);
    if isCropXY
        nn = 8; % size param
        xx1 = 451; xx2 = xx1+2^nn-1; % crop 
        yy1 = 451; yy2 = yy1+2^nn-1; % crop
        save('cropCoor','xx1','xx2','yy1','yy2');
    else
        xx1 = 1; xx2 = En1; % crop
        yy1 = 1; yy2 = Boy1; % crop
    end
    [Boy1,En1]=size(imgFrst(yy1:yy2,xx1:xx2));
    imageInfo=imfinfo(fname);
    Frames=length(imageInfo);
    %Frames = 10; % # of frames to be processed

    coeffMat = dir('coeff*.mat');
    isCoeffLoaded = 0; Coeff = 0;
    if numel(coeffMat)==1
        load(coeffMat.name)
        display('loading coefficient')
        isCoeffLoaded = 1;
    elseif numel(coeffMat) > 1
        display('ERROR : more than one coeff.mat files');
        return;
    end
    assignFileNames

    %% FIND X & Y
    if ~exist(posDataFileNm) 
        nG = 5; ng = (nG-1) /2;
        if isPixelsNonIsotropic
            nG2 = 7; ng2 = (nG2-1) /2;
            x = -ng2:ng2; y = x;
            [X Y] = meshgrid(x,y);
            xyRat = dx/dy; % x-y pixel size ratios
            x2 = [-1:1]*xyRat; y2 = -ng:ng;
            [X2 Y2] = meshgrid(x2,y2);
            gaus=fspecial('gaussian', nG2, 1);
            gaus = interp2(X,Y,gaus,X2,Y2);
            %gaus = gaus(:,2:4);
        else
            gaus=fspecial('gaussian', nG, 1);
        end
        lap=[-1,-1,-1;-1,8,-1;-1,-1,-1];
        
        if ~isCoeffLoaded
            tit = 'determine the coefficient for background filtering';
            hFig=figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256));
            set(gcf,'units','normalized','outerposition',[0 0 1 1])
            kk =0.1;
            axe=axes('Parent',hFig,'DataAspectRatio',[1 1 1],'Position',[kk kk 1-2*kk 1-2*kk ],'Visible','off','outerposition',[0 0 1 1]);

            % TileParameters
            mag = 2; % magnification
            temp=get(0,'ScreenSize');
            screenX = temp(3); screenY = temp(4);
            nX = round((screenX + 200)/En1/mag);
            nY = round((screenY - 0)/Boy1/mag);
            nFrTile = nX*nY; % number of frames to be tiled
            nRow = nY;  
            
            for j=1:nFrTile
                temp=imread(fname,j);
                PREdata(:,:,j) = temp(yy1:yy2,xx1:xx2);
                PREdataFilt = imfilter(double(PREdata),gaus,'symmetric');
                PREdataFilt = imfilter(PREdataFilt,lap,'symmetric');
            end
            PREdataTiled = tileFrames(PREdata ,nRow);
            PREdataFiltTiled = tileFrames(PREdataFilt ,nRow);

            SHOW=PREdataTiled;
            imgMx = max(PREdataFiltTiled(:));
            Coeff = imgMx/5;

            hTextCoeff1 = uicontrol('style','text','BackgroundColor',[1 1 1],'String',sprintf('Coefficient: %i',round(Coeff)),'Position',[20 150 160 15]);
            hTextCoeff2 = uicontrol('style','text','BackgroundColor',[1 1 1],'String',{'select the coefficient'; 'use keys to change threshold value';'w: ++0.1 e:--0.1';'s: ++0.01 d:--0.01';'press q to continue'},'Position',[20 20 160 115]);
            q=0; 
            slideVal = double(double(Coeff)/double(imgMx));
            while q == 0
                colormap(pink);
                imagesc(PREdataTiled);axis image;
                PREdataFiltTiledDiv=PREdataFiltTiled/Coeff;
                BINAR=im2bw(PREdataFiltTiledDiv,1);
                DINAR=uint16(BINAR).*PREdataTiled;
                BWSHOW=imregionalmax(DINAR, 8);
                [y,x,v]=find(BWSHOW==1);
                figure(1);
                %imagesc(BWSHOW); % filtered spots
                hold on; scatter(x,y,22,'o');
                hold off;
                set(hTextCoeff1,'String',sprintf('Coefficient: %i',round(Coeff)));
                axis image;
                btn = 0;
                while btn == 0
                    btn = waitforbuttonpress;
                    k = get(hFig,'CurrentCharacter');
                end
                upd= 0;
                switch lower(k)
                    case 's'
                        slideVal = slideVal - 0.01; upd = 1;
                    case 'd'
                        slideVal = slideVal + 0.01; upd = 1;
                    case 'w'
                        slideVal = slideVal - 0.1; upd = 1;
                    case 'e'
                        slideVal = slideVal + 0.1; upd = 1;
                    case 'q'
                        q = 1;
                end
                 if upd == 1
                     Coeff = slideVal*imgMx;
                 end
            end
            save('Coeff.mat','Coeff');
            assignFileNames
        end
        close;
        
        for j=1:Frames
            img =imread(fname,j);
            J(:,:,j) = img(yy1:yy2,xx1:xx2);
            img = double(img(yy1:yy2,xx1:xx2));
            
            JF(:,:,j) = imfilter(img,gaus,'symmetric');
            JF(:,:,j) = imfilter(JF(:,:,j),lap,'symmetric');
        end
        J=uint16(J);
        h = waitbar(0,'3D localization...');
        %% find 3D intensity
        sp = 0;
        nSp = 0;
        Coeff = double(Coeff);
        for k=1:Frames
            if exist('BREAK.MAT')
                if exist('breakVal.mat')
                    clear
                    load breakVal;
                    delete('BREAK.MAT','breakVal')
                else
                    save('breakVal');
                    return;
                end
            end
            IMG = J(:,:,k);
            % find peaks
            dataFilt = JF(:,:,k); 
            dataFiltDiv = dataFilt/Coeff;
            bin = im2bw(dataFiltDiv,1);
            if isempty(find( bin == 1)),continue;end
            din = uint16(bin).*IMG;
            BW = imregionalmax(din, 8);
            [B,L] = bwboundaries(BW,'noholes');

            q=0; 
            %DEFINE Size
            %[En,Boy]=size(IMG);
            nSp = nSp + length(B);
            for m=1:length(B) % for each spot
                c=cell2mat(B(m));
                %csize=(max(c(:,1))-min(c(:,1)))*(max(c(:,2))-min(c(:,2)));
                q=q+1;sp = sp + 1;
                Py=uint16(mean(c(:,1)));
                Px=uint16(mean(c(:,2)));
                Px0 = double(Px);
                Py0 = double(Py);
                
                % adjust the big window center position for the spots at
                % the edges
                PxBW = Px;
                PyBW = Py;
                     
                if (Px-(BigWindowSizeX+1)/2)<1
                    PxBW=(BigWindowSizeX+1)/2;
                end
                if (Py-(BigWindowSizeY+1)/2)<1
                    PyBW=(BigWindowSizeY+1)/2;
                end
                if (Px+(BigWindowSizeX+1)/2)>En1
                    PxBW=En1-(BigWindowSizeX+1)/2;
                end
                if (Py+(BigWindowSizeY+1)/2)>Boy1
                    PyBW=Boy1-(BigWindowSizeY+1)/2;
                end
                dPy = double(PyBW)-double(Py);
                dPx = double(PxBW)-double(Px);
                Px = PxBW; Py = PyBW;
                %DEFINE Window
                yy = Py-(WindowSizeY+1)/2;
                y1 = yy + 1; y2 = yy + WindowSizeY;
                xx = Px-(WindowSizeX+1)/2;
                x1 = xx + 1; x2 = xx + WindowSizeX;
                Window=IMG(y1:y2,x1:x2);

                %DEFINE Big Window
                yy = PyBW-(BigWindowSizeY+1)/2;
                y1 = yy + 1; y2 = yy + BigWindowSizeY;
                xx = PxBW-(BigWindowSizeX+1)/2;
                x1 = xx + 1; x2 = xx + BigWindowSizeX;
                BigWindow=IMG(y1:y2,x1:x2);

                %background
                BACKmean=[min(mean(BigWindow,1)),min(mean(BigWindow,2))];
                BACK(q,k)=min(BACKmean);

                %FIND Total Intensity
                INT(q,k)=sum(sum(Window))-BACK(q,k)*(WindowSizeX*WindowSizeY);

                %FIND Intensity Center
                TopX=0;
                TopY=0;
                TopColum=0;
                TopRow=0;
                WSumX=0;
                WSumY=0;

                for j=1:WindowSizeX
                   TopX(j)=sum(Window(:,j));
                end
                TopX=TopX-min(TopX);
                TopRow=sum(TopX);

                for j=1:WindowSizeX
                    WSumX=WSumX+j*TopX(j);
                end

                for i=1:WindowSizeY
                   TopY(i)=sum(Window(i,:));
                end
                TopY=TopY-min(TopY);
                TopColum=sum(TopY);

                for i=1:WindowSizeY
                    WSumY=WSumY+i*TopY(i);
                end

                Xc(k)=WSumX/TopRow;
                Yc(k)=WSumY/TopColum;
                if isnan(Xc(k)), Xc(k)=(WindowSizeX+1)/2; end;
                if isnan(Yc(k)), Yc(k)=(WindowSizeY+1)/2; end;


                PXc=uint8(Xc(k));
                PYc=uint8(Yc(k));

                %center of intensity
                X_=double(Px)+Xc(k)-double((WindowSizeX+1)/2);
                Y_=double(Py)+Yc(k)-double((WindowSizeY+1)/2);
        
                X(q,k)=X_-dPx;
                Y(q,k)=Y_-dPy;
                %Inten(Py,Px)=INT(k);
                         
                if StackNum > 1
                    %% FIND Z NOW
                    for l=1:StackNum
                        DigitDiff=floor(log10(Frames))-floor(log10(k));   

                        if k == 1
                        DigitDiff=floor(log10(Frames));
                        end

                        if DigitDiff == 0   
                        Stack=imread(['stack_' int2str(k) '.tif'],l);
                        end

                        if DigitDiff == 1   
                        Stack=imread(['stack_0' int2str(k) '.tif'],l);
                        end

                        if DigitDiff == 2   
                        Stack=imread(['stack_00' int2str(k) '.tif'],l);
                        end

                        if DigitDiff == 3   
                        Stack=imread(['stack_000' int2str(k) '.tif'],l);
                        end

                        Stack=double(Stack);
                        Stack=Stack-min(min(Stack));

                        ZWindow(:,:,l)=Stack(Py-(WindowSizeY+1)/2+1:Py+(WindowSizeY+1)/2-1,Px-(WindowSizeX+1)/2+1:Px+(WindowSizeX+1)/2-1);

                    end
                    TopZr=squeeze( sum(sum(ZWindow,2),1)-sum(min(ZWindow)*WindowSizeX,2));
                    % older ver. : TopZr(l)=sum(sum(ZWindow))-sum(min(ZWindow))*WindowSize; %  data range in the XY plane
                    zr = 1:StackNum;

                    z=zr;
                    TopZ=TopZr';
                    TopZ=TopZ-min(TopZ); % remove background
                    % fz=fit(z,TopZ,'gauss1');
                    % coeffvalues(fz);
                    % Z(q,k)=ans(2);
                    Z(q,k)=StackNum-sum(z.*TopZ)/sum(TopZ);

                    zPosZero = 0;
                    if TopZ==0
                        %fprintf('z position is zero check the code');
                        zPosZero = zPosZero + 1;
                        Z(q,k)=StackNum/2;
                        %return;
                    end
                else
                    Z = zeros(size(X));
                end
            end
            waitbar(k/Frames)
            if ~rem((k-1),100) 
                time100frame = 0;
                if k > 1
                    time100frame = toc; 
                    disp(sprintf('time for 100 frames : %.02f\n',time100frame));
                end;
                tic;
                
            end
        end
        disp(sprintf('rate of gaussian fit success : %.02f',numel(find(X>0))/nSp));
        %stackWrite(NBINs*3000,'NBINs.tif')
        clear JF J ;
        close(h);
        save(xyzDataFileNm,'X','Y','Z','INT'); 
        clear X Y Z INT spotWin NBINs Zstack maxProj;
        clear NBIN NDIN BW B L IMG NBINsum temp imgFrst Frames Frames2 ZWindow Stack
        save(posDataFileNm); % save data  

        
        %% move files
        folderNm= sprintf('%s-coeff%i',label,round(Coeff));
        if exist(folderNm)
            fn = 1;
            while exist(sprintf('%s-coeff%i_%i',label,round(Coeff),fn))
                fn = fn + 1;
            end
            folderNm = sprintf('%s-coeff%i_%i',label,round(Coeff),fn);
            mkdir(folderNm);
        else
            mkdir(folderNm);
        end

        % move files
        isCopy = 1;
        if isCopy
            copyfile('coeff.mat',strcat(folderNm,'\',nmCoeff));
            copyfile(fname,strcat(folderNm,'\',fname));            
            if exist('cropCoor.mat'), movefile('cropCoor.mat',strcat(folderNm,'\cropCoor.mat')); end;
            movefile(posDataFileNm,strcat(folderNm,'\',posDataFileNm));
            movefile(xyzDataFileNm,strcat(folderNm,'\',xyzDataFileNm));
            if exist('fname.mat')
                copyfile('fname.mat',strcat(folderNm,'\fname.mat'));         
            else
                copyfile('inputInfo_wesScope.mat',strcat(folderNm,'\inputInfo_wesScope.mat'));         
            end
            cd(folderNm)
        end
        %return;
    end

    %% TIME To FIND OUT THE TRACES
    if ~exist(traceDataFileNm0)
        load(posDataFileNm);
        load(xyzDataFileNm);

        isCropFrames = 0;
        if isCropFrames
            f1=1;f2=300;
            X = X(:,f1:f2);
            Y = Y(:,f1:f2);
            Z = Z(:,f1:f2);
        end
        Xilk=X;
        Yilk=Y;
        Zilk=Z;
        [Boy,Frames]=size(X);

        h = waitbar(0,'Finding the traces...');
        p=0;
        f=0;
        tic
        isDebug = 0;
        for k=1:Frames-1 % number of frames
            if isDebug, disp(sprintf('=frame#:%i\n',m)); end
            for m=1:Boy % number of spots
                if X(m,k)>0
                if X(m,k)<En1
                    tracex=zeros(1,Frames);
                    tracey=zeros(1,Frames);
                    tracez=zeros(1,Frames);
                    traceint=zeros(1,Frames);

                    dif=Inf(Boy,Frames-k+1);
                    difbin=zeros(Boy,Frames-k+1);            
                    AslX=X(m,k);       % last X value in the tracking
                    AslY=Y(m,k);       % last Y value in the tracking
                    AslZ=Z(m,k);
                    if isDebug, disp(sprintf('==spot#:%i\n',m)); end
                    quit = 0;
                    ll=1;
                    for l=1:Frames-k+1 % later frames                        
                        if isDebug, disp(sprintf('===check frame#:%i\n',l)); end
                        if quit, 
                            break; 
                        end
                        for n=1:Boy % all spots
                            dif(n,l)=sqrt((AslX-X(n,k+l-1))^2*dx^2 + (AslY-Y(n,k+l-1))^2*dy^2+(AslZ-Z(n,k+l-1))^2*dz^2);        
                        end
                        [v n2] = min(dif(:,l));
                        BOY = Boy;
                        BOY = 1;
                        for n=1:BOY % all spots
                            n=n2;
                            if l-ll>sptReAppearTime % if the next frame where a spots re-appears in sptJmp distance is 4 frames apart ignores it.
                                %if sum(sum(difbin(:,l-2:l-1))) == 0, 
                                    quit =1;
                                    break, 
                                %end
                            end
                            
                            if isDebug, disp(sprintf('===check spot#:%i\n',n)); end
                            if sqrt((AslX-X(n2,k+l-1))^2+(AslY-Y(n2,k+l-1))^2+[(AslZ-Z(n2,k+l-1))*PlaneDist/PixelSize]^2) < sptJmp;
                                if dif(n,l)==min(dif(:,l))
                                    if n ~= n2
                                        disp(sprintf('n:%i, n2"%i\n',n,n2));
                                        %disp('strange')
                                        break
                                    end         
                                    ll=l;
                                    difbin(n,l)=1;
                                    AslX=Xilk(n,k+l-1);
                                    AslY=Yilk(n,k+l-1);
                                    AslZ=Zilk(n,k+l-1);
                                    break
                                end
                            end
                        end % spots
                    end % frames

                    for n=1:Boy
                    for l=1:Frames-k+1
                        if difbin(n,l)==1;
                            tracex(k+l-1)=Xilk(n,k+l-1);
                            tracey(k+l-1)=Yilk(n,k+l-1);
                            tracez(k+l-1)=Zilk(n,k+l-1);
                            traceint(k+l-1)=INT(n,k+l-1);
                            X(n,k+l-1)=Inf;
                            Y(n,k+l-1)=Inf;
                            Z(n,k+l-1)=Inf;
                        end
                    end % frames
                    end % spots

                    num=numel(find(tracex>0));
                    if num>=minTraceLength % if # of data points larger than minTraceLength, than saves as a trace
                         pos=find(tracex>0);
                         ilk=zeros(1,num+1);
                         son=zeros(1,num+1);
                         ilk(1:1:num)=pos(1:1:num);
                         son(2:1:num+1)=pos(1:1:num);
                         fark=ilk-son;
                         %if numel(find(fark==1))>2 % # of consecutive data points
                        p=p+1;
                        TraceX(p,:)= sparse(tracex);
                        TraceY(p,:)=sparse(tracey);
                        TraceZ(p,:)=sparse(tracez);
                        TraceINT(p,:)=sparse(traceint);   
                        

                         %end
                    end
                end % x<EN
                end % x>0
            end
           waitbar(k / Frames)
           if exist('_stopRunning-ON')
               break
           end
        end
        TrackTime = toc;
        save('TrackTime','TrackTime');
        close(h)
        %return;
        save(traceDataFileNm0,'TraceX','TraceY','TraceZ','TraceINT'); 
        save(cfgTraceFileNm,'sptJmp','sptReAppearTime','minTraceLength');
        delete(traceDataFileNm)
    end
    
    if ~exist(traceDataFileNm) % speed and trace combination
        %% COMBINE TRACES
        load(traceDataFileNm0)
        [m n]=size(TraceX);
        h = waitbar(0,'Combining the traces...');
        disp('Combining the traces...');

        [Boy2,Frames]=size(TraceX);
        tic;
        com=0;
        for i=1:Boy2 % each trace
        %     TraceDif1=0;TraceDif2=0;TraceDif3=0;
            if 0 && i <= Boy2-1 && isCombTraces  % combine traces
                LastElement=max(find(TraceX(i,:)>0));
                LastBefore=LastElement-1;
                for j=i+1:Boy2 
                    FirstElement=min(find(TraceX(j,:)>0));
                    FirstAfter=FirstElement+1;
                    if FirstElement-LastElement>-1
                        if FirstElement-LastElement<3
                            TraceDif1=sqrt([TraceX(i,LastElement)-TraceX(j,FirstElement)]^2+[TraceY(i,LastElement)-TraceY(j,FirstElement)]^2+[(TraceZ(i,LastElement)-TraceZ(j,FirstElement))*PlaneDist/PixelSize]^2);
                            TraceDif2=sqrt([2*TraceX(i,LastElement)-TraceX(i,LastBefore)-TraceX(j,FirstElement)]^2+[2*TraceY(i,LastElement)-TraceY(i,LastBefore)-TraceY(j,FirstElement)]^2+[(2*TraceZ(i,LastElement)-TraceZ(i,LastBefore)-TraceZ(j,FirstElement))*PlaneDist/PixelSize]^2);
                            TraceDif3=sqrt([TraceX(i,LastElement)-2*TraceX(j,FirstElement)+TraceX(j,FirstAfter)]^2+[TraceY(i,LastElement)-2*TraceY(j,FirstElement)+TraceY(j,FirstAfter)]^2+[(TraceZ(i,LastElement)-2*TraceZ(j,FirstElement)+TraceZ(j,FirstAfter))*PlaneDist/PixelSize]^2);
                            TraceDif=[TraceDif1,TraceDif2,TraceDif3];
                            %TraceDif=sqrt((TraceX(i,LastElement)-TraceX(j,FirstElement))^2+(TraceY(i,LastElement)-TraceY(j,FirstElement))^2);
                            if min(TraceDif)<10
                                for k=1:Frames
                                    if [TraceINT(i,k)+TraceINT(j,k)] > 0
                                    TraceX(i,k)= [TraceX(i,k)*TraceINT(i,k)+TraceX(j,k)*TraceINT(j,k)]/[TraceINT(i,k)+TraceINT(j,k)];
                                    TraceY(i,k)= [TraceY(i,k)*TraceINT(i,k)+TraceY(j,k)*TraceINT(j,k)]/[TraceINT(i,k)+TraceINT(j,k)];
                                    TraceZ(i,k)= [TraceZ(i,k)*TraceINT(i,k)+TraceZ(j,k)*TraceINT(j,k)]/[TraceINT(i,k)+TraceINT(j,k)];
                                    TraceINT(i,k)= [TraceINT(i,k)*TraceINT(i,k)+TraceINT(j,k)*TraceINT(j,k)]/[TraceINT(i,k)+TraceINT(j,k)];
                        %             TnaceX(i,k)=TraceX(i,k); TnaceY(i,k)=TraceY(i,k); TnaceX(j,k)=TraceX(j,k); TraceY(j,k)=TraceY(j,k);
                                    end
                                end
                                TraceX(j,:)=NaN; TraceY(j,:)=NaN; TraceZ(j,:)=NaN; TraceINT(j,:)=NaN;
                                LastElement=max(find(TraceX(i,:)>0));
                                LastBefore=LastElement-1;
                                com=com+1;
                            end
                        end
                    end
                end
            end % i <= Boy2-1
                %Tcomp0 = toc;
            %% speed
            % fill gaps (missing frames in jumpy traces)
            tracex = TraceX(i,:);
            tracey = TraceY(i,:);
            tracez = TraceZ(i,:);
            if isnan(tracex(1)), continue; end;
            
            ind = find(tracex~=0);
            frst = min(ind); last = max(ind);
            jmp = abs((tracex(frst:last)>0)-1);
            if sum(jmp > 0)
                bnd = bwboundaries(jmp,'noholes');
                for j = 1:numel(bnd) % for each jump
                    temp = bnd{j}+frst-1; % jump boundaries
                    jb = temp(:,2);
                    mx= max(jb); mn=min(jb);
                    jL = mx-mn+2; % length

                    jSx = tracex(mx+1)-tracex(mn-1); % size
                    jSy = tracey(mx+1)-tracey(mn-1); % size
                    jSz = tracez(mx+1)-tracez(mn-1); % size
                    jsX = jSx/jL;% step
                    jsY = jSy/jL;% step
                    jsZ = jSz/jL;% step
                    jVx = tracex(mn-1)+jsX*(1:jL-1);
                    jVy = tracey(mn-1)+jsY*(1:jL-1);
                    jVz = tracez(mn-1)+jsZ*(1:jL-1);
                    tracex(mn:mx)=jVx;
                    tracey(mn:mx)=jVy;
                    tracez(mn:mx)=jVz;
                    
                end
                TraceX(i,:) = tracex;
                TraceY(i,:) = tracey;
                TraceZ(i,:) = tracez;
            end                       
            % speed
            xDiff = [tracex(2:end)-tracex(1:end-1) 0]*PixelSize;
            yDiff = [tracey(2:end)-tracey(1:end-1) 0]*PixelSize;
            zDiff = [tracez(2:end)-tracez(1:end-1) 0]*PlaneDist;
            traceSpeed = sqrt(xDiff.^2 + yDiff.^2 + zDiff.^2);   
            
            ind = find(tracex~=0);
            frst = min(ind); last = max(ind);
            if frst ~= 1, traceSpeed(frst-1) = 0; end
            if last ~= Frames ,traceSpeed(last) = 0; end;
            TraceSpeed(i,:) = traceSpeed;
                
        end
        waitbar(i / Boy2)
        close(h)
        % 'traceData-coeff%d.mat'
        Tcomp = toc;
        save(traceDataFileNm,'TraceX','TraceY','TraceZ','TraceINT','TraceSpeed','Tcomp')
        
    end


%% PLOT TRACES
    disp('preparing plot data...')
    if isCombTraces
        load(traceDataFileNm);
    else
        load(traceDataFileNm0);
    end
    figure
    hist(TraceSpeed((TraceSpeed(:)>0) .* (TraceSpeed(:)<20)>0),1000);
    
    %TraceSpeed = round(TraceSpeed*100)/100;
    %TraceSpeed(TraceSpeed~=3.17) = 0.1;
    
    % clear traces on the edge
    TraceX = full(TraceX);
    TraceY = full(TraceY);
    TraceZ = full(TraceZ);
    TraceSpeed = full(TraceSpeed);
    TraceX(TraceX==0) = nan;
    TraceY(TraceY==0) = nan;
    TraceZ(TraceZ==0) = nan;
    TraceSpeed(TraceSpeed==0) = nan;
    
    mxX = max(TraceX,[],2);
    mnX = min(TraceX,[],2);
    mxY = max(TraceY,[],2);
    mnY = min(TraceY,[],2);
    discardTraces = find((mnX<=5) + (mxX>En1-4) + (mnY<=5) + (mnY>Boy1-4));
    
    for i = 1:numel(discardTraces)
        TraceX=TraceX([1:discardTraces(i)-1 discardTraces(i)+1:end],:);
        TraceY=TraceY([1:discardTraces(i)-1 discardTraces(i)+1:end],:);
        TraceZ=TraceZ([1:discardTraces(i)-1 discardTraces(i)+1:end],:);
        TraceSpeed=TraceSpeed([1:discardTraces(i)-1 discardTraces(i)+1:end],:);
        discardTraces = discardTraces - 1;
    end
    
    % Tracespeed mean
    if 0
        for i = 1:size(TraceSpeed,1)
            TraceSpeed(i,TraceSpeed(i,:)>0)= mean(TraceSpeed(i,TraceSpeed(i,:)>0));
        end
    end
    
    nAv=3;
    for i = 1:size(TraceSpeed,1)
        TraceSpeed(i,TraceSpeed(i,:)>0)= conv(TraceSpeed(i,TraceSpeed(i,:)>0),ones(nAv,1),'same');
    end
    
    
    % color data
    CplotVecN = size(TraceX,1); % # of traces
    Nframe = size(TraceX,2); % # of frames
    CMinZ = 1; % oclor map in Z
    CMinSpeed =0; % color in speed
    %speedMax = 3.5;
    if CMinZ 
        Zmin = min(TraceZ(TraceZ~=0));
        Zmax = max(TraceZ(:));
        Zrange = Zmax - Zmin;
        Cplot = round(63*(TraceZ-Zmin) / Zrange)+1;
        CM = colormap('jet');
    elseif CMinSpeed
        Zmin = min(TraceSpeed(TraceSpeed~=0));
        Zmax = max(TraceSpeed(:));
        Zrange = Zmax - Zmin;
        Cplot = round(63*(TraceSpeed-Zmin) / Zrange)+1;
        CM = colormap('jet');
    else % color by trace
        CplotVec = 1: CplotVecN;
        CplotVec = rem(CplotVec,64)+1;
        Cplot = repmat(CplotVec',1,Nframe);
        CM = colormap('lines');
    end
    Frames = size(TraceX,2);
    En = Frames;
    [Boy2]=size(TraceX,1);
    
    nonZeroIx = find(TraceY>0);
    dispPx = {'*',2}; % plot3k
    minimum = 1;
    maximum = max(Cplot(:));
    
%% PLOT
    NAN = find(isnan(TraceX));
    TraceX(NAN)=0;
    TraceY(NAN)=0;
    TraceZ(NAN)=0;

    % parameters:
    load(posDataFileNm,'sptReAppearTime'); %(frames) use the value from tracker function generating trace values
   
    imgZFout = 'TraceImage';
    magImg = 6;

    img2D = imread(fname,1); 
    img2D = img2D(yy1:yy2,xx1:xx2);
    imSz = size(img2D');
    TraceY = imSz(2)-TraceY+1;
    
    figSz(1) = imSz(1)*magImg;
    figSz(2) = imSz(2)*magImg;
    colormap('gray');

    zMax = max(TraceZ(:));
    if zMax == 0, zMax = 1; TraceZ = TraceZ+1; end;

    % find the frames where the traces disappear
    hat = ones(1,sptReAppearTime); 
    dspTrcFrm = zeros(size(TraceX,1),1);
    for i = 1:size(TraceX,1) % all traces
        trace = full(TraceX(i,:));
        trFrm = conv(trace,hat,'same'); % active frames of the trace
        trFrm(1:find(trFrm>0,1)) = 1;
        lastFrm = find(trFrm==0,1); 
        if ~isempty(lastFrm)
            dspTrcFrm(i) = lastFrm;
        end
    end

    imgZFout = 'TraceImage';
    hQ = 0; %hImg = image; 
    lastX = TraceX(:,1); 
    lastY = TraceY(:,1);
    lastZ = TraceZ(:,1);
    tit = 'image';
    m = imSz(1); n = imSz(2);
    pos=get(0,'ScreenSize');
    pos=pos(3:4) - [m n-35];
    fig=figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos/2 m n]);
    axe=axes('Parent',fig,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren','XLim',0.5+[0 m],'YLim',0.5+[0 n]);

    for ixFrm = 1:En-1
        img2D = imread(fname,ixFrm+1); 
        img2D = flipud(img2D(yy1:yy2,xx1:xx2));
        %delete(hImg);

        hImg = imagesc(img2D,'Parent',axe); %axis image; 
        %set(gca,'position',[0 0 1 1]); set(gcf,'position',[0 0 figSz(1) figSz(2)]);
        ixTrc = find(TraceX(:,ixFrm+1)>0);
        currX = TraceX(:,ixFrm+1);
        currY = TraceY(:,ixFrm+1);
        currZ = TraceZ(:,ixFrm+1);
        uistack(hImg,'bottom'); hold on
        ixTrc2 = find(lastX(ixTrc)>0);
        ix = ixTrc(ixTrc2); % indices of traces tracked in the current frame    
        dspTrcIx = find(ixFrm+1 == dspTrcFrm); % index for dissappearing traces
        showTrace =1;
        if showTrace
            hQdel = hQ(dspTrcIx,:); % handles for discont. traces
            delete(hQdel(find(hQdel~=0)));    % remove the traces of the discontinued traces
        end
        
        for i = 1:round(length(ix)) % for each trace 
            iL = ix(i);  % index for each line
            hQ(iL,ixFrm)=quiver(lastX(iL),lastY(iL),currX(iL)-lastX(iL),currY(iL)-lastY(iL),'Color',CM(Cplot(iL,ixFrm),:));
            if ~showTrace
                adjust_quiver_arrowhead_size(hQ(iL,ixFrm-1),5)
            end
        end
        % print images
        imgFig = getframe(gcf);
        imgOut = imgFig.cdata;
        figPos = get(gcf,'Position');
        %imgOut = imgOut((figPos(4)-m)/2+1:(figPos(4)+m)/2,(figPos(3)-n)/2+1:(figPos(3)+n)/2,:);
        imgOut = imgOut((figPos(4)-n)/2+1:(figPos(4)+n)/2,(figPos(3)-m)/2+1:(figPos(3)+m)/2,:);
        fOut = sprintf('%s.tif',imgZFout);
        if ixFrm == 1
            if exist([imgZFout '.tif'])
                delete([imgZFout '.tif']);
            end
            imwrite(imgOut,fOut,'Compression', 'none') 
        else
            imwrite(imgOut,fOut,'WriteMode','append','Compression', 'none') 
        end

        delete(hImg)
        %export_fig(imgZFout,'-append');
        fprintf('frame %i/%i \n',ixFrm,En);
        lastX = currX; 
        lastY = currY;
        lastZ = currZ;
    end
    plotColorBar(TraceSpeed/acqRate_sec,6,'speed [\mum/sec]');
    set(gcf,'PaperPositionMode','auto');
    colorbarPrint = sprintf('%s_colorbar.tiff',fOut);
    print(colorbarPrint,'-dtiff','-r80'); 


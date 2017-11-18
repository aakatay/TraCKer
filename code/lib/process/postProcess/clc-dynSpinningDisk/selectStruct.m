% selects structures from two channel movie (dyn-clc)
% run in the cell folder
% dyn sub folder with dyn movie
clear

pdiffThresh = 300;


xytLim = [136 333 156 377];
xytLim = [134 250 144 264];
xytLim = [278 196 292 204];

xytLim = [179 181 190 194 54]; % ix1
%xytLim = [320 161 330 174 58]; %ix2
%xytLim = [310 96 320 114 67]; % ix3
%xytLim = [235 150 250 170 45]; % ix4
%xytLim = [225 140 240 160 58]; % ix5
% xytLim = [255 100 270 110 30]; % ix6 2nd best

% 260,104,32 cell5
% 232,153,61
% 243 164 46
% 313 102 68
% 324 169 59
% 184, 185 , 54
% filter out some (isDisp =1)
id = [160 165 167 168 187 193 197]; % CELL1
id = []; % cell2
%id = [228]; % cell5

%skipFrames0 = [25.5]; % cell2
%skipFrames0 = [36.5]; % cell2
%skipFrames = unique([ceil(skipFrames0)-[0:3] floor(skipFrames0)+[0:3]]);
%skipFrames = [1:11 126:258]; % cell1
skipFrames = [];

isSelectStructs = 0; % set isReAdjustProfLines=1
isPits = 0; % process pits

isDispCombined = 0; % multiple cells
isDisp =0;

isReAdjustProfLines = 1; 
skipGoodOnes = 0;
isIxRedo = 0;

tDiffSingleFrame = 0;

if isPits
    if exist('profData.mat'), warning('rename old profData file'); pause; end;
    fid = fopen('pits.txt');
    i=1;
    while 1
        pitsIx_ = fgets(fid);
        if pitsIx_ == -1, break; end;
        pitsIx(i) = str2num(pitsIx_);
        i=i+1;
    end
    fclose(fid);
end

pl0 = 41; % profile length
profC0=nan(1,pl0);
profC2=nan(1,pl0);
profClate=nan(1,pl0); 
profD0=nan(1,pl0);
profD2=nan(1,pl0);
profDlate=nan(1,pl0); 
if isIxRedo, load('ixRedo'); end


if ~isDispCombined
    fnTandemData = rdir('tandemTrack.mat');
    load(fnTandemData(1).name); % 'xc','yc','xd','yd','fc','fd'
    xxc0 = xc;
    yyc0 = yc;
    xxc = round(xc);
    yyc = round(yc);
    xxd0 = xd;
    yyd0 = yd;
    xxd = round(xd);
    yyd = round(yd);
    
    

    fnacqDIR = rdir('clc\acq_*');
    fnacq = fnacqDIR.name;
    fnacq2DIR = rdir('dyn\acq_*');
    fnacq2 = fnacq2DIR.name;
    dirName = pwd; % directory name
    dS = strfind(dirName,'\'); % dirSlash
    dc=2;
    last = numel(dirName);
    if dc~=0, last=dS(end-dc+1)-1; end
    cellName = dirName(dS(end-dc)+1:last);
    dateName = dirName(dS(end-1-dc)+1:dS(end-1-dc)+6);
    cellLabel = [dateName '-' cellName];
    fnameClcDyn = sprintf('%s-ClcDynOverlay.tif',cellLabel);
    fnameDynTime = sprintf('%s-DynFrameNo.tif',cellLabel);
    profPlotFN  = sprintf('%s-profilePlots.tif',cellLabel);
    isFramesDefined = 0;
    %dt = 1; % # of frames after the peak to com
    dt = 3;

    if exist(profPlotFN), delete(profPlotFN); end;

    Cimg = imread(fnacq,1); % clathrin channel
    


    imginfo = imfinfo(fnacq);
    frm2=numel(imginfo);

    fs = 9; % frame size
    fs2 = 4; % frame size
    nd = numel(xxc) ; % number of colocalized diff peaks
    
    for i = 1:frm2 % read acquisition images
        A_ = imread(fnacq,i); % clathrin channel
        A(:,:,i) = A_;
        D_ = imread(fnacq2,i); % dynamin channel
        D(:,:,i) = D_;
    end
end

if ~isDisp

    figure(101);clf
    figure(100);clf
    figure(97); clf
%fclose('all');delete('ROI.txt')
%close all

    %fndyn = rdir('dyn\_**\*traceData0*.mat');
    %load(fndyn.name); % trInf


    if exist('frames.txt')
        isFramesDefined = 1;
        fid = fopen('frames.txt');
        inp = fgetl(fid);
        sc = strfind(inp,':'); % pos. semicolon
        frm1 = str2num(inp(1:sc-1));
        frm2 = str2num(inp(sc+1:end));
        ndel = sum(trInf(:,1)>frm2);
        trInf(end-ndel+1:end,:)=[];
        fclose(fid);
    end



    Dimg = zeros(size(Cimg));
    Dtimg = zeros(size(Cimg));
    Dtimg(sub2ind(size(Cimg),yyd,xxd))=fc; % peak times

    % intensity of dynamin
    Dimg(sub2ind(size(Cimg),yyd,xxd))=1;  
    figure(1); imagesc(Dimg);axis image


    isLoadProfData = 0;
    if exist('profData.mat'),isLoadProfData=1; load profData; end;

    for i = 1:nd % find structs
        
        
        if fc(i)<xytLim(5), continue; end;
        
        if isPits,  if ~ismember(i,pitsIx), continue; end; end
        isContinue=0;
        %if i<189, continue; end;
        
        if isIxRedo, if ~ismember(i,ixRedo), continue; end; end
        dispImages;
        
        if isContinue, continue; end;
        
        drawProfile; % p0_1 p0_2 --> prof0, prof2
        if isContinue, continue; end;

        % adjust profile
        figure(100); title(sprintf('struct#:%i (frm#:%i)',i,tc));
        him = imline(gca,[e1x e2x],[e1y e2y]);
        hold on ; hs=scatter(e1x,e1y,'r','o','filled'); 
        scatter(xd0,yd0,'r','.'); 
        scatter(xc0,yc0,'c','.')
        hold off;
        set(gca,'Xlim',xl)
        set(gca,'Ylim',yl)
        P0=getPosition(him);
        while 1 % adjust line position
            figure(100)
            waitforbuttonpress
            pause(0.4)
            P1=getPosition(him);
            Pc = sum((P0-P1)~=0,2); % detect change
            Pce = find(Pc>0); % index of the end changed position
            if isempty(Pce) % no changes
                kp=get(gcf,'CurrentCharacter')+1; 
                if kp == 28 % esc key (skip)
                    e1x=nan; e1y=nan; e2x=nan; e2y=nan;
                    profC0(i,:)=nan;
                    profC2(i,:)=nan;
                    profClate(i,:)=nan;
                    profD0(i,:)=nan;
                    profD2(i,:)=nan;
                    profDlate(i,:)=nan;
                    break
                elseif kp == 14 % enter key (save)
                    break
                elseif kp == 114 % q key (save and quit)
                    if ~isReAdjustProfLines
                        profC0 = profC0(1:end-1,:);
                        profC2 = profC2(1:end-1,:);
                        profClate = profClate(1:end-1,:);
                        profD0 = profD0(1:end-1,:);
                        profD2 = profD2(1:end-1,:);
                        profDlate = profDlate(1:end-1,:);
                    end
                    save('profData','exy','profC0','profC2','profClate','profD0','profD2','profDlate'); % 0: on the edge, nan: not good shape
                    return;
                end
                continue;
            elseif isempty(find(Pc==0)) % vector moved
                P1 = P0; % ignore changes
                Pce = 2;
            end            
            cx = P1(Pce,1);
            cy = P1(Pce,2);


            ex = (xd0-cx)*2+cx; % endpoints
            ey = (yd0-cy)*2+cy;
            e1x=ex; e2x=cx; e1y=ey; e2y=cy;
            Pce2 = 2-Pce+1;
            P2(Pce2,:) = [e1x e1y];
            P2(Pce,:) = [e2x e2y];
            delete(him); if exist('hs'), delete(hs); end
            him = imline(gca,[P2(1) P2(2)],[P2(3) P2(4)]);
            hold on ; hs=scatter(P2(1),P2(3),'r','o','filled'); hold off;
            P0=getPosition(him);
            pause(0.1)
            drawProfile; % e1x e2x,e1y e2y --> prof0, prof2
        end
        exy(i,:) = [e1x e1y e2x e2y];

    end
    exy(sum(exy,2)==0,:)=nan;
    save('profData','exy','profC0','profC2','profClate','profD0','profD2','profDlate'); % 0: on the edge, nan: not good shape
end
exy_ = [];
profC0_=[];
profC2_=[];
profClate_=[];
profD0_=[];
profD2_=[];
profDlate_=[];
nProf_ = [];
if isDispCombined % multiple cells
    profDIR  = rdir('*\profDataSel.mat');
    for i = 1:numel(profDIR)
        load(profDIR(i).name);
        nProf_(i) = size(exy,1);
        exy_ = [exy_;exy];
        profC0_ = [profC0_;profC0];
        profC2_ = [profC2_;profC2];
        profClate_ = [profClate_;profClate];
        profD0_ = [profD0_;profD0];
        profD2_ = [profD2_; profD2];
        profDlate_ = [profDlate_; profDlate];
    end
    exy = exy_;
    profC0=profC0_;
    profC2=profC2_;
    profClate=profClate_;
    profD0=profD0_;
    profD2=profD2_;
    profDlate=profDlate_;
    if isempty(exy), error('trying to save empty array'); end
    nProf = nProf_;
    save('profData','nProf','exy','profC0','profC2','profClate','profD0','profD2','profDlate'); % 0: on the edge, nan: not good shape
else
    load('profData');
end


%% display results
figure(201);
delete(201)
figure(201); maximize
%ixp = find(sum(profC0>0,2)>0); % defined profiles
ixp = find(sum(exy>0,2)>0); % defined profiles

% filter out some
ixp(find(ismember(ixp,id)))=[];


np = numel(ixp);
profFilt=nan(size(profC0));
profFilt(ixp,:)=1;

plen = sum(profC0_>=0,2);

hs = uicontrol('Style', 'slider','Max',np,'Position',[20 20 600 20],'SliderStep', [1/(np),1/(np)]);
ht = uicontrol('Style','text','Position',[20 60 90 20]);
pC0 = profC0.*profFilt;
pC2 = profC2.*profFilt;
pCl = profClate.*profFilt;
pD0 = profD0.*profFilt;
pD2 = profD2.*profFilt;
pDl = profDlate.*profFilt;
ixs0 = 1;
ixs = 0;
pl0 = 41;
isQuit = 0;
while 1
    while ixs0 == ixs
        ixs = round(get(hs,'Value')); % selected index
        kp=get(gcf,'CurrentCharacter')+1; 
        if kp == 28 % esc key (quit)
            isQuit = 1;
            break;
        end
        pause(0.2)
    end
    if isQuit, break; end;
    
    ixs0 = ixs;
    
    clf(201)
    
    isPlotAll=0;
    hs = uicontrol('Style', 'slider','Max',np,'Position',[20 20 600 20],'Value',ixs,'SliderStep', [1/(np),1/(np)]);
    ht = uicontrol('Style','text','Position',[20 60 90 20]);
    if ixs %ixs >0
        ixdel = find(ixs~=ixp);
        set(ht,'String',int2str(ixp(ixs)));
    else 
        ixdel = find(~ismember(1:44,ixp));
        isPlotAll=1;
    end
        
%     pC0(ixde
%%
% $x^2+e^{\pi i}$ l,:)=[];
%     pC2(ixdel,:)=[];
%     pD0(ixdel,:)=[];
%     pD2(ixdel,:)=[];
    
    vx = [1:pl0]-21;
    subplot(4,2,1); hp(1,:)=plot(vx,pC0'); title('clathrin, C0, (t=t0)'); grid minor;
    subplot(4,2,3); hp(2,:)=plot(vx,pC2'); title('clathrin, C1, (t=t0+1)'); grid minor;
    subplot(4,2,5); hp(3,:)=plot(vx,(pC0-pC2)'); title('clathrin diff, C0-C1 (t0+1 - t0)'); grid minor;
    subplot(4,2,2); hp(4,:)=plot(vx,pD0'); title('dynamin, D0, (t=t0)'); grid minor;
    subplot(4,2,4); hp(5,:)=plot(vx,pD2'); title('dynamin, D1, (t=t0+1)'); grid minor;
    subplot(4,2,6); hp(6,:)=plot(vx,(pD0-pD2)'); title('dynamin diff, D0-D1 (t0+1 - t0)'); grid minor;

    CL = jet(np);
    for i = 1:np % set color gradings
        LW = 0.5;
        if i == ixs, LW = 2; end
        set(hp(1,ixp(i)),'Color',CL(i,:),'LineWidth',LW);
        set(hp(2,ixp(i)),'Color',CL(i,:),'LineWidth',LW);
        set(hp(3,ixp(i)),'Color',CL(i,:),'LineWidth',LW);
        set(hp(4,ixp(i)),'Color',CL(i,:),'LineWidth',LW);
        set(hp(5,ixp(i)),'Color',CL(i,:),'LineWidth',LW);
        set(hp(6,ixp(i)),'Color',CL(i,:),'LineWidth',LW);
        
        if i == ixs
            uistack(hp(1,ixp(i)), 'top');
            uistack(hp(2,ixp(i)), 'top');
            uistack(hp(3,ixp(i)), 'top');
            uistack(hp(4,ixp(i)), 'top');
            uistack(hp(5,ixp(i)), 'top');
            uistack(hp(6,ixp(i)), 'top');
        end
    end
    
    if isPlotAll  
        mpC0 = mean(pC0(ixp,:),1);
        vx(isnan(mpC0))=[];
        mpC0(isnan(mpC0))=[];
        mpC2 = mean(pC2(ixp,:),1);mpC2(isnan(mpC2))=[];
        mpD0 = mean(pD0(ixp,:),1);mpD0(isnan(mpD0))=[];
        mpD2 = mean(pD2(ixp,:),1);mpD2(isnan(mpD2))=[];
        mpC02 = (mpC0-mpC2);
        mpD02 = (mpD0-mpD2);
        pcComp = [mpC0; mpC2; mpC02*3]';
        pdComp = [mpD0; mpD2; mpD02*2]';
        pcComp2 = mean(pCl(ixp,:),1)';pcComp2(isnan(pcComp2))=[];
        pdComp2 = mean(pDl(ixp,:),1)';pdComp2(isnan(pdComp2))=[];
        title7 = sprintf('clathrin mean N=%i',np);
        title8 = sprintf('dynamin mean N=%i',np);
    else
        pcComp = [pC0(ixp(ixs),:); pC2(ixp(ixs),:); pC0(ixp(ixs),:)-pC2(ixp(ixs),:)]';
        pdComp = [pD0(ixp(ixs),:); pD2(ixp(ixs),:); pD0(ixp(ixs),:)-pD2(ixp(ixs),:)]';
        pcComp2 = pCl(ixp(ixs),:)';
        pdComp2 = pDl(ixp(ixs),:)';

        i = ixp(ixs);
        isReAdjustProfLines = 1;
        isContinue=0;
        isLoadProfData = 1;
        dispImages;
        drawProfile;
        title7 = sprintf('(x,y,frm)=(%.02f,%.02f,%i)',xc(i),yc(i),fc(i));
        title8 = sprintf('(x,y,frm)=(%.02f,%.02f,%i)',xd(i),yd(i),fd(i));
    end
    figure(201)
    subplot(4,2,7); 
    %plotyy(vx,pcComp,vx,pcComp2); title(sprintf('clathrin mean N=%i',np)); grid minor;
    plot(vx,pcComp); title(title7); grid minor;
    legend('C0','C1','C0-C1(x3)')
    subplot(4,2,8); plot(vx,pdComp); title(title8); grid minor;
    legend('D0','D1','D0-D1(x2)','Dlate')
    
    if isDispCombined
        % save to txt
        xlsFN = 'profileCompare.xls';

        % convert to cells
        v1s = mpC0.*0+1; % vector of ones
        mpC0c = mat2cell(mpC0,1,v1s);
        mpC2c = mat2cell(mpC2,1,v1s);
        mpC02c = mat2cell(mpC02,1,v1s);
        pcComp2c = mat2cell(pcComp2',1,v1s);
        mpD0c = mat2cell(mpD0,1,v1s);
        mpD2c = mat2cell(mpD2,1,v1s);
        mpD02c = mat2cell(mpD02,1,v1s);
        pdComp2c = mat2cell(pdComp2',1,v1s);

        %write to XLS
        xlsData1 = [[{'C0'} mpC0c]; [{'C1'} mpC2c];[{'Cdiff'} mpC02c];[{'Clate'} pcComp2c]];
        xlsData2 = [[{'D0'} mpD0c]; [{'D1'} mpD2c];[{'Ddiff'} mpD02c];[{'Dlate'} pdComp2c]];
        xlsData = [xlsData1;xlsData2];
        xlswrite(xlsFN,xlsData);
    end
    
    
    ndp = numel(mpC0);
    
    
    %pause(0.3)
end
delete(hs); delete(ht);
profCompFN = 'profileCompare.tif';
imgFig = getframe(gcf);
dataImg = imgFig.cdata; 
imwrite(uint16(dataImg),profCompFN); % structMap (selected structures)


exy = exy(ixp,:);
profC0 = profC0(ixp,:);
profC2 = profC2(ixp,:);
profClate = profClate(ixp,:);
profD0 = profD0(ixp,:);
profD2 = profD2(ixp,:);
profDlate = profDlate(ixp,:);  
save('profDataSel','nProf','exy','profC0','profC2','profClate','profD0','profD2','profDlate'); % 0: on the edge, nan: not good shape



return;


































%%
% ---------
    if 0
        imgFig = getframe(gcf);
        dataImg = imgFig.cdata; 
        imwrite(uint16(dataImg),profPlotFN,'WriteMode','append','Compression', 'none'); % 

        % draw on dyn-clc image
        figure(102)
        mx = max([A0(:)' A2(:)'] );
        clear imgCol
        imgCol(:,:,1)=A2/mx; % red
        imgCol(:,:,2)=A0/mx; % green
        imgCol(:,:,3)=A2.*0;
        imshow(imgCol, 'InitialMagnification','fit');  axis image; % maximize;
        hold on
        line([e1x e2x],[e1y e2y],'Color','r');
        scatter(xd,yd,'r','*')
        scatter(cx,cy,'b','.')
        hold off

        % dynamin intensity 3x3 average vs time
        figure(103) % dynamin intensity vs time
        plot(tvec,Dtprf)
        hold on;
        line([t-t0+1 t-t0+1],[0 max(Dtprf)],'Color','r');
        hold off;
        grid
    end
% ---------

% overlay image
mn = min(Cimg(:));
mx = max(Cimg(:));
CM = parula(mx);
CM(end+1,:)= [1 0 0 ]; % red
CDimg = Cimg;
CDimg(sub2ind(size(Cimg),yy,xx)) = mx+1;
figure(5)
imagesc(CDimg); axis image
colormap(CM);
set(gcf,'units','pixels','Position',[120,120,size(Cimg,2),size(Cimg,1)]); 
set(gca,'units','pixels','Position',[0,0,size(Cimg,2),size(Cimg,1)]);
imgFig = getframe(gcf);
dataImg = imgFig.cdata; 
imwrite(uint16(dataImg),fnameClcDyn); % structMap

% ---------



if ~isFramesDefined % select time window
    figure(2); imagesc(Dtimg);axis image
    CM = jet(max(tfc));
    colormap(CM);
    set(gcf,'units','pixels','Position',[120,120,size(Cimg,2),size(Cimg,1)]); 
    set(gca,'units','pixels','Position',[0,0,size(Cimg,2),size(Cimg,1)]);
    imgFig = getframe(gcf);
    dataImg = imgFig.cdata; 
    imwrite(uint16(dataImg),fnameDynTime); % structMap
    figure(3)
    hist(tfc) 
    disp('set last frame')
    return;
end
% --------- 

 div = 1;
    WindowSize = 5; 
    BigWindowSize=WindowSize;
    En1 = size(D,2);
    Boy1 = size(D,1);
    clear Dtprf;
    for j=1:size(Dt,3)
        Px = round(xdg);
        Py = round(ydg);
        IMG = D(:,:,tvec(j));
        centOfMassLoc; % Px Py --> INT_
        Dtprf(j) = INT_;
    end
% ---------
   % auto find threshold for boundary detection
        figure(99) % clathrin pixel intensity histogram
        nbins=round(numel(A0)/4/3);
        hist(A0(:),nbins)
        [bins,val]=hist(A0(:),nbins);
        threshold = val(find(max(bins)==bins)+2);
        threshold = threshold(end);
        A0f = A0; 
        A0f(A0<threshold)=0; % filtered
        %A0f = zeros(size(A0f));        A0f(3:4,3:4)=1;
        [Bt,MASK]=bwboundaries(A0f,'noholes');
        [max_size, max_index] = max(cellfun('size', Bt, 1));
        Bt = Bt{max_index};
        MASK(MASK~=max_index)=0;
        A0f=A0f.*MASK;
        szA = size(A0f);
        sy = szA(1); sx = szA(2);
        Ax = repmat([1:sx],sy,1);
        Ay = repmat([1:sy]',1,sx);
        Axw = A0f.*Ax;
        cx = mean(Axw(:))/mean(A0f(:));
        Ayw = A0f.*Ay;
        cy = mean(Ayw(:))/mean(A0f(:));
        % display boundaries
        figure(100)
        imagesc(A0); axis image;
        hold on;
        bp = plot(Bt(:,2), Bt(:,1), 'r', 'LineWidth', 2);
        scatter(cx,cy,'.','r')
        scatter(xd,yd,'r','*')
        hold off;
        
        
        
        % set boundaries(select threshold)
        isbreak=0;isEdit=0;
        pixdel = ones(size(A0)); % delete pixels
        while (1) % set coeff 
            while ~isEdit
                figure(100)
                kp=get(gcf,'CurrentCharacter')+1;
                if kp ==31 % up key
                    threshold = threshold + 25;
                    break
                elseif kp ==32 % down key
                    threshold = threshold - 25;
                    break
                elseif kp == 28 % esc key = > edit mode
                    isEdit = 1;
                    break
                elseif kp == 14 % enter key  = > done
                    isbreak=1;
                    break
                end
                pause(0.1)
            end
            if 0&isEdit % delete pixels mode
                figure(100)
                [gx,gy] = ginput(1);
                pixdel(round(gy),round(gx))=0;
                if isempty(gx)
                    kp=get(gcf,'CurrentCharacter')+1;
                    if kp == 14 % enter key  = > done
                        isEdit = 0;
                    end
                end
            end
            A0f = A0; 
            A0f(A0<threshold)=0;% filtered
            A0f = A0f.*pixdel;
            [Bt,MASK]=bwboundaries(A0f,'noholes');
            [max_size, max_index] = max(cellfun('size', Bt, 1));
            Bt = Bt{max_index};
            MASK(MASK~=max_index)=0;
            A0f=A0f.*MASK;
            
            delete(100)
            figure(100)
            imagesc(A0.*pixdel); axis image;
            title(sprintf('coeff:%0.1f',threshold));
            hold on;
            plot(Bt(:,2), Bt(:,1), 'r', 'LineWidth', 2)
            Axw = A0f.*Ax;
            cx = mean(Axw(:))/mean(A0f(:));
            Ayw = A0f.*Ay;
            cy = mean(Ayw(:))/mean(A0f(:));
            hold on 
            scatter(cx,cy,'.','r')
            scatter(xd,yd,'r','*')
            hold off;
            if isbreak, 
                stArr{i} = pixdel.*(A0>=threshold);
                break; 
            end
        end
        
        
        
% ------

        if 0 
            figure(21); % 2-color image
            mags=2;
            mx = max([A0(:)' A2(:)'] );
            clear imgCol
            imgCol(:,:,1)=A2/mx; % red
            imgCol(:,:,2)=A0/mx; % green
            imgCol(:,:,3)=A2.*0;
            %imgCol(fs+1,fs+1,:)=[0 0 1];
            imsz=size(A0);
            imsz = imsz*mags;
            imshow(imgCol);  axis image; % maximize;
            %imshow(repelem(imgCol,mags,mags,1));  axis image; % maximize;
            hold on;
            scatter(xc0,yc0,'.','b');
            hold off;
            set(gca,'Units','pixels'); 
            set(gca,'Position',[0 0 imsz])
            set(gcf,'Units','pixels'); 
            set(gcf,'Position',[200 200 imsz])

            figure(22)
            imagesc(A0)        
            hold on;
            scatter(xd,yd,'.','b');
            hold off;
            set(gca,'Units','pixels'); 
            set(gca,'Position',[0 0 imsz])
            set(gcf,'Units','pixels'); 
            set(gcf,'Position',[200 200 imsz])
        end
% ------

        
        %%
        plen = sum(profC0>=0,2);
        ixRedo = find(((plen<=11) .* (plen~=0)))
        ixRedo = find(((plen<=11) .* ~isnan(mean(exy,2))))
        
        %ixRedo=[6,10,13,18,19];
        save('ixRedo','ixRedo')
        
        %%
% ------
clear
load profData
i = 1:216;
exy(i,:) = nan;
profC0(i,:) = nan;
profC2(i,:) = nan;
profClate(i,:) = nan;
profD0(i,:) = nan;
profD2(i,:) = nan;
profDlate(i,:) = nan;    
save('profData','exy','profC0','profC2','profClate','profD0','profD2','profDlate'); % 0: on the edge, nan: not good shape


%% ----- save all open figs
figHandles = get(0,'Children');
for i = 1:numel(figHandles)
    figure(i);
    imgFig = getframe(gcf);
    dataImg = imgFig.cdata; 
    FN=sprintf('%i.tif',i);
    imwrite(uint16(dataImg),FN); % structMap (selected structures)
    
end
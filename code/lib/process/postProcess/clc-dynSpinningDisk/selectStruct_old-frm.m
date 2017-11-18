% selects structures from two channel movie (dyn-clc)
% run in the cell folder
% dyn sub folder with dyn movie

%fclose('all');delete('ROI.txt')
%close all
figure(101);clf
figure(100);clf
figure(97); clf
clear
fnacqDIR = rdir('acq1_*');
fnacq = fnacqDIR.name;
fnacq2DIR = rdir('dyn\acq_*');
fnacq2 = fnacq2DIR.name;
dirName = pwd; % directory name
dS = strfind(dirName,'\'); % dirSlash
cellName = dirName(dS(end)+1:end);
dateName = dirName(dS(end-1)+1:dS(end-1)+6);
cellLabel = [dateName '-' cellName];
fnameClcDyn = sprintf('%s-ClcDynOverlay.tif',cellLabel);
fnameDynTime = sprintf('%s-DynFrameNo.tif',cellLabel);
profPlotFN  = sprintf('%s-profilePlots.tif',cellLabel);
isFramesDefined = 0;
%dt = 1; % # of frames after the peak to com
dt = 3;

if exist(profPlotFN), delete(profPlotFN); end;

mag = 1; % magnification
%magOld =2;
Cimg = imread(fnacq,1); % clathrin channel
Cimg = repelem(Cimg,mag,mag);

fndyn = rdir('dyn\_**\*traceData0*.mat');
load(fndyn.name); % trInf

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


tf = trInf(:,1)+trInf(:,2)-1; % frame number
xx0 = trInf(:,4)*mag;
yy0 = trInf(:,5)*mag;
xx = round(xx0);
yy = round(yy0);
int = trInf(:,6);
Dimg = zeros(size(Cimg));
Dtimg = zeros(size(Cimg));
Dtimg(sub2ind(size(Cimg),yy,xx))=tf;

% intensity of dynamin
Dimg(sub2ind(size(Cimg),yy,xx))=int;  
figure(1); imagesc(Dimg);axis image

if ~isFramesDefined % select time window
    figure(2); imagesc(Dtimg);axis image
    CM = jet(max(tf));
    colormap(CM);
    set(gcf,'units','pixels','Position',[120,120,size(Cimg,2),size(Cimg,1)]); 
    set(gca,'units','pixels','Position',[0,0,size(Cimg,2),size(Cimg,1)]);
    imgFig = getframe(gcf);
    dataImg = imgFig.cdata; 
    imwrite(uint16(dataImg),fnameDynTime); % structMap
    figure(3)
    hist(tf) 
    disp('set last frame')
    return;
end

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

nd = size(trInf,1); % number of dyn peaks
for i = 1:frm2 % find structs around dyn peaks
    A_ = imread(fnacq,i); % clathrin channel
    A(:,:,i) = repelem(A_,mag,mag);
    D_ = imread(fnacq2,i); % dynamin channel
    d(:,:,i) = D_;
    D(:,:,i) = repelem(D_,mag,mag);
    
end

fs = 7; % frame size
if ~exist('ROI.txt') % make selections of structures
    fid = fopen('ROI.txt','wt');
    
    for i = 1:nd % find structs
        t=tf(i);
        x = xx(i);
        y = yy(i);
        if round(x)<2 | round(y)<2 | round(x)>size(Cimg,2)-2 | round(y)>size(Cimg,1)-2
            continue;
        end
if y<40, continue; end; % false recs
        xf = xx0(i); % double precision
        yf = yy0(i);
        x0 = x-fs; if x0<1, x0=1; end
        x1 = x+fs; if x1>size(Cimg,2), x1=size(Cimg,2); end;
        y0 = y-fs; if y0<1, y0=1; end
        y1 = y+fs; if y1>size(Cimg,1), y1=size(Cimg,1); end;
        t0 = t-2; if t0<1, t0=1; end;
        t1 = t+1; if t1>frm2, t1=frm2; end;
        t2 = t+dt; if t2>frm2, t2=frm2; end;
        
        mags=20;
        
        xd = xf-x0+0.5;
        yd = yf-y0+0.5;
        xd = xd;
        yd = yd;

        A0 = mean(A(y0:y1,x0:x1,t0:t),3);
        A2 = mean(A(y0:y1,x0:x1,t1:t2),3);
        mx = max([A0(:)' A2(:)'] );
        clear imgCol
        imgCol(:,:,1)=A2/mx; % red
        imgCol(:,:,2)=A0/mx; % green
        imgCol(:,:,3)=A2.*0;
        %imgCol(fs+1,fs+1,:)=[0 0 1];
        imsz=size(A0);
        imsz = imsz*mags;
        
        figure(21)
        imshow(imgCol);  axis image; % maximize;
        %imshow(repelem(imgCol,mags,mags,1));  axis image; % maximize;
        hold on;
        scatter(xd,yd,'.','b');
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
        
        sel(i)=1;
        [~, r]=imcrop(gcf); % rectangle pos
        if isempty(r)
            sel(i)=0;
            continue
        end
        r(1:2)= floor(r(1:2)/2)*2+1;
        r(3:4)= ceil(r(3:4)/2)*2;
        r=r+[x0-1 y0-1 0 0];
        fprintf(fid,'%iX%iY%ix%ifrm%idyn%.02fx%.02fy\n',[r t xf yf]);
    end
    save('selectStruct','sel');
    fclose(fid)
    return;
else % process structures
    fid = fopen('ROI.txt','r');
    fid2 = fopen('ROI2.txt','wt');
    pl0 = 41; % profile length
    prof0 = zeros(1,pl0);
    prof2 = zeros(1,pl0);
    i= 1;
     
    isLoadProfData = 0;
    if exist('profData.mat'),isLoadProfData=1; end;
    
    while 1 % read all coors of structs
        %close all;
        
        
        coors = fgetl(fid);
        if sum(coors == -1) || isempty(coors), break, end % eof
        if isempty(str2num(coors(1))) % not a stable structure
            coors = coors(2:end);
            isNonStruct = 1;
            continue;
        end        
        
        
        % skip till
        if i < 0,
            i = i + 1;
            continue;
        end
        
        
        % read ROI.txt for struct coors and pos
        a=sscanf(coors,'%dX%dY%dx%dfrm%ddyn%fx%fy');
    a([1:4])=round(a([1:4]));
     a(6:7) = a(6:7);
        disp([coors '_Ix' int2str(i)])
        xdg = a(6); % global coors
        ydg = a(7);
        x0 = a(1); y0 = a(2);  x1 = x0+a(3)-1; y1 = y0 + a(4)-1; 
        t = a(5); xd = xdg-x0+1; yd = ydg-y0+1;
        yd0 = round(round(ydg)/2)*2-mag-1; yd1 = round(round(ydg)/2)*2+mag; 
        xd0 = round(round(xdg)/2)*2-mag-1; xd1 = round(round(xdg)/2)*2+mag;
        
        t0 = t-8; if t0<1, t0=1; end;
        t0_2 = t-1; if t0_2<1, t0_2=1; end;
        t1 = t+1; if t1>frm2, t1=frm2; end;
        t2_1 = t+1; if t2_1>frm2, t2_1=frm2; end;
        t2_2 = t+2; if t2_2>frm2, t2_2=frm2; end;
        t2_3 = t+3; if t2_3>frm2, t2_3=frm2; end;
        
        % clathrin intensity averages at diff. times around dyn peak
        A0 = mean(A(y0:y1,x0:x1,t0:t),3);
        A0_2 = mean(A(y0:y1,x0:x1,t0_2:t),3);
        A2 = mean(A(y0:y1,x0:x1,t1:t2_1),3);
        A3 = mean(A(y0:y1,x0:x1,t2_2:t2_2),3);
        A4 = mean(A(y0:y1,x0:x1,t2_3:t2_3),3);
        D0 = mean(D(y0:y1,x0:x1,t0:t),3);
        D0_2 = mean(D(y0:y1,x0:x1,t0_2:t),3);
        D2 = mean(D(y0:y1,x0:x1,t1:t2_1),3);
        D3 = mean(D(y0:y1,x0:x1,t2_2:t2_2),3);
        D4 = mean(D(y0:y1,x0:x1,t2_3:t2_3),3);
        tvec = t0:t2_2+5;
        Dt = D(yd0:yd1,xd0:xd1,tvec);
        nm=numel(Dt(:,:,1));
        
        
        div = 1;
        WindowSize = 5; 
        BigWindowSize=WindowSize;
        En1 = size(d,2);
        Boy1 = size(d,1);
        clear Dtprf;
        for j=1:size(Dt,3)
            Px = round(xdg);
            Py = round(ydg);
            IMG = D(:,:,tvec(j));
            centOfMassLoc; % Px Py --> INT_
            Dtprf(j) = INT_;
        end
        %Dtprf = mean(reshape(Dtconv,nm,numel(Dt)/nm),1);
        
        %% images
        figure(97) % clathrin images before and after dyn peak
        subplot(4,3,1)
        imagesc(A0); axis image; title('A0'); 
        subplot(4,3,2)
        imagesc(A2); axis image; title('A2'); 
        subplot(4,3,3)
        imagesc(A3); axis image; title('A3'); 
        subplot(4,3,4)
        imagesc(D0); axis image; title('D0'); 
        subplot(4,3,5)
        imagesc(D2); axis image; title('D2'); 
        subplot(4,3,6)
        imagesc(D3); axis image; title('D3'); 
        
        subplot(4,3,7)
        imagesc(A0); axis image; title('A0'); 
        subplot(4,3,8)
        imagesc(uint16(A0-A2)); axis image; title('Adiff2'); 
        subplot(4,3,9)
        imagesc(uint16(A2-A3)); axis image; title('Adiff3'); 
        subplot(4,3,10)
        imagesc(D0); axis image; title('D0'); 
        subplot(4,3,11)
        imagesc(uint16(D0-D2)); axis image; title('Ddiff2'); 
        subplot(4,3,12)
        imagesc(uint16(D2-D3)); axis image; title('Ddiff3'); 
        
        %%
        figure(100)
        imagesc(A0); axis image;
        
        if isLoadProfData
            ex=exy(i,1);
            ey=exy(i,2);
            cx=exy(i,3);
            cx=exy(i,4);
        else
        
            cx = xd - 1; cy = yd;
            ex = (xd-cx)*2+cx; % endpoints
            ey = (yd-cy)*2+cy;
            e1x=cx; e2x=ex; e1y=cy; e2y=ey;
        end
        % xy line profile
        [px1,py1,p0_1] = improfile(A0_2,[e1x xd],[e1y yd]);
        [px2,py2,p0_2] = improfile(A0_2,[xd e2x],[yd e2y]);
        
        drawProfile;
        figure(100)
        him = imline(gca,[e1x e2x],[e1y e2y]);
        P0=getPosition(him);
        while 1
            hold on; scatter(xd,yd,'r','*'); hold off;
            P1=getPosition(him);
            Pc = sum((P0-P1)~=0,2); % change
            Pce = find(Pc>0); % index of the end changed position
            if isempty(Pce) % no changes
                pause(0.4)
                kp=get(gcf,'CurrentCharacter')+1; 
                if kp == 28 % esc key
                    break
                end
                continue;
            end            
            cx = P1(Pce,1);
            cy = P1(Pce,2);
                
            ex = (xd-cx)*2+cx; % endpoints
            ey = (yd-cy)*2+cy;
            e1x=cx; e2x=ex; e1y=cy; e2y=ey;
            Pce2 = 2-Pce+1;
            P2(Pce,:) = [e1x e1y];
            P2(Pce2,:) = [e2x e2y];
            delete(him);
            him = imline(gca,[P2(1) P2(2)],[P2(3) P2(4)]);
            P0=getPosition(him);
            pause(0.1)
            drawProfile; % --> prof0, prof2
        end
        exy(i,:) = [e1x e1y e2x e2y];
        

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
        	
        i = i + 1;
        delete(100)
    end
    save('profData','exy','prof0','prof2');
    fclose(fid)
    mx = max(prof0,[],2);
    prof0n = prof0./repmat(mx,[1,pl0]);
    prof2n = prof2./repmat(mx,[1,pl0]);
    plot(prof0n','b')
    hold on 
    plot(prof2n','r')
    hold off
    line([ceil(pl0/2) ceil(pl0/2)],[0,1.5])
    
end

return;
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

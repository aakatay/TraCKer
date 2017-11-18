function [f1,f2,binImgLastTime,binImgCr] = cropTimeWinCore(binFN,f1,f2,cP)
% crops recruimentTrack data in time
% called by cropTimeWin & find structs
% cP : cropPos
    binTime = 5; % [sec]

    %CM = jet(nf);
    if ~exist(binFN), binFN = ['../' binFN]; end
    nf = numel(imfinfo(binFN));
    figure(13);
    set(gcf,'name',sprintf('time crop, frames:%i-%i',1,nf))

    for i = 1:nf
        binImg(:,:,i) = imread(binFN,i);
        binImgTime(:,:,i) = im2bw(binImg(:,:,i),0)*i;
    end
    binImg = binImg(cP(2):cP(4),cP(1):cP(3),:);
    binImgTime = binImgTime(cP(2):cP(4),cP(1):cP(3),:);
    
    

        % recruitmentRate plot
        figure(14); 
        set(gcf,'name','set time window');
        sx = size(binImg,2);
        sy = size(binImg,1);
        sz = nf;
        recruitmentRate = sum(reshape(double(binImg),sx*sy,sz),1)/5; % per sec
        % calculate half bleach time
        mxrecRate=max(recruitmentRate);
        recruitmentRateSmth = smooth(recruitmentRate);
        hlfRecRate=mxrecRate/2;
        ix1 = find(recruitmentRateSmth< hlfRecRate,1);
        v1=recruitmentRateSmth(ix1-1);
        v2=recruitmentRateSmth(ix1);
        hlfBlchTime = (ix1+(v2-hlfRecRate)/(v1-v2))*binTime;
        %plot
        xTime = (1:sz); % [sec]
        plot(xTime,recruitmentRate); 
        ax = gca;
        ax.XTick = [(0:3:sz)];
        grid minor;
        ylabel('# of recruitment incidences/sec')
        %xlabel('time [sec]');
        xlabel('frames');
        title(sprintf('recruitmentRate[incidence/sec]. half bleach time: %.02fsec',hlfBlchTime))

    Rc = binImg(:,:,1); 
    set(0,'units','pixels');szScr = get(0,'screensize');
    szScr = szScr(3:4);
    szRc = size(Rc);
    mag=round(min(szScr./szRc/6));
    if mag<1, mag = 1; end;
    colormap('jet')
    cset = 0;
    while (1) % set time window (enter to quit)
        if f2 > nf, f2 = nf; end;
        binImgCr = sum(binImg(:,:,f1:f2),3);
        binImgLastTime = max(binImgTime(:,:,f1:f2),[],3);
        figure(13)


        imagesc(binImgLastTime); 
        set(gcf,'name',sprintf('time crop, frames:%i-%i',f1,f2))
        axis image;
        if cset
            set(gca,'Xlim',xl)
            set(gca,'Ylim',yl)
        else
            set(gcf,'units','pixels','Position',[50,szScr(2)-size(Rc,1)*mag-100,size(Rc,2)*mag,size(Rc,1)*mag]); 
            set(gca,'units','pixels','Position',[0,0,size(Rc,2)*mag,size(Rc,1)*mag]);
        end;

        figure(14);
        p1 = (f1);
        p2 = (f2);
        hl1 = line([p1 p1],[0 mxrecRate]);
        hl2 = line([p2 p2],[0 mxrecRate]);
        [X,Y] = ginput(1);
        cset = 1;
        delete(hl1);delete(hl2);
        if isempty(X)
            break;
        end
        if abs(X-f1) > abs(X-f2)
            f2 = round(X);
        else
            f1 = round(X);
        end
        figure(13) 
        xl = get(gca,'Xlim');
        yl = get(gca,'Ylim');
    end
end
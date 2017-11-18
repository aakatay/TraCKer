% generates data for the scatter plot of prebleach intensity vs recruitment
% numbers
clear all;
isdbg = 1;
dirCell = rdir('**\cell*');
nc = numel(dirCell);

cellDir = 1;
if isempty(dirCell)
    cellDir =0;
    nc=1;
end
PWD = pwd;
for j = 1:nc
    close all;

    if cellDir, cd(dirCell(j).name); end
    %%

        dirName = pwd; % directory name
        dS = strfind(dirName,'\'); % dirSlash
        cellName = dirName(dS(end)+1:end);
        dateName = dirName(dS(end-1)+1:dS(end-1)+6);
        cellLabel = [dateName '-' cellName];
    %%

    avgFN = rdir('Avg_lap_*.tif');
    recFN = rdir('binImgRcrtSum_time*.tif');

    stFN = rdir('recComp\structMapUpd*.mat');
    pxStatFN = 'pxStat.mat';
    pxStatScatFN = sprintf('pxStatScat_%s.tif',cellLabel);

    load(stFN.name)

    A = imread(avgFN.name);
    A = repelem(A,4,4);
    R = imread(recFN.name);

    stn = size(LtNew,3);
    as = [];
    rs = [];
    Ac = {};
    Rc = {};
    for i = 1:stn
        lt = LtNew(:,:,i);
        lt(lt==0)=nan;
        a = lt.*double(A);
        r = lt.*double(R);
        if isdbg
            figure(2)
            imagesc(a)
            figure(3)
            imagesc(r)
            figure(4)
            imagesc(lt)
            cc=3;
        end
        
        a=a(:);
        a(isnan(a))=[];
        r=r(:);
        r(isnan(r))=[];
        a=a/(max(a));
        r=r/(max(r));
        Ac{i} = a;
        Rc{i} = r;
        as = [as; a];
        rs = [rs; r];
    end


    rv = 1:max(rs)+1;
    rsv = unique(rs);
    amn = [];
    for i = 1:numel(rsv)
        amn(i) = mean(as(rs==rsv(i)));
    end

    figure(1)
    hold on;
    dx=linspace(-0.4,0.4,stn)*0;
    for i = 1:stn
        scatter(Rc{i}+dx(i),Ac{i},[],ones(numel(Ac{i}),1)*i,'.');
    end
    plot(rsv,amn,'k')
    hold off;
    colormap('jet')

    set(gcf,'color','w');

    XLIM = [-0.2 max(rs)+0.2];
    xlim(XLIM)
    ylabel('intensity')
    xlabel('recruitment number')
    maximize
    title(cellLabel)

    imgFig = getframe(gcf); 
    imgOut = imgFig.cdata;
    imwrite(imgOut,pxStatScatFN,'Compression', 'none') 

    save(pxStatFN,'Rc','Ac')
    cd(PWD)
end
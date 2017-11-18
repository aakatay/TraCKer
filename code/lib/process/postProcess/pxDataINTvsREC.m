% generates data for the scatter plot of prebleach intensity vs recruitment
% numbers
clear all;
isCompileALL = 0;
nhc = 20; % number of hist counts
isdbg = 1;

if ~isCompileALL
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
        pxStatFN = 'pxStat2.mat';
        pxStatScatFN = sprintf('pxStatScat2_%s.tif',cellLabel);

        load(stFN.name)

        A = imread(avgFN.name);
        A = repelem(A,4,4);
        R = imread(recFN.name);

        stn = size(LtNew,3);
        stn = 4;
        as = [];
        rs = [];
        Ac = {};
        Rc = {};
        for i = stn:stn % each structure
            lt = LtNew(:,:,i);
            lt(lt==0)=nan;
            a = lt.*double(A);
            r = lt.*double(R);
            
             a(:,all(isnan(a)))=[];a(all(isnan(a),2),:)=[];
             r(:,all(isnan(r)))=[];r(all(isnan(r),2),:)=[];
            if isdbg
                figure(2)
                imagesc(a)
                figure(3)
                imagesc(r)
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
            pause;
        end


        [N,edges]=histcounts(as,nhc);
        rmn = [];
        asv = [];
        for i = 1:numel(edges)-1
            rmn(i) = mean(rs(  find((edges(i)<as) .* (as<edges(i+1))) ));
            asv(i) = (edges(i)+edges(i+1) )/2;
        end

        figure(1)
        hold on;
        dx=linspace(-0.4,0.4,stn)*0;
        for i = 1:stn % each structure
            scatter(Ac{i}+dx(i),Rc{i},[],ones(numel(Ac{i}),1)*i,'.');
        end
        plot(asv,rmn,'k')
        hold off;
        colormap('jet')

        set(gcf,'color','w');
        xlabel('intensity')
        ylabel('recruitment number')
        maximize
        title(cellLabel)

        imgFig = getframe(gcf); 
        imgOut = imgFig.cdata;
        imwrite(imgOut,pxStatScatFN,'Compression', 'none') 

        save(pxStatFN,'Rc','Ac')
        cd(PWD)
        pause
    end
else

    %% compile all
    close all;

    pxStatALLFN = 'pxStat2ALL.mat';
    pxStatScatALLFN = 'pxStatScat2ALL.tif';
    pxStatScatALLstackedFN = 'pxStatScat2ALLstacked.tif';

    pd = rdir('*\**\pxStat2.mat'); % px data from all cells
    imgFN = rdir('*\**\pxStatScat2_*.tif'); % px data from all cells

    AcAll = {};
    RcAll = {};
    for i = 1:numel(pd)
        load(pd(i).name);
        AcAll = [AcAll Ac];
        RcAll = [RcAll Rc];
        imgOut=imread(imgFN(i).name);
        imwrite(imgOut,pxStatScatALLstackedFN,'WriteMode','append','Compression', 'none') 
    end

    AcAll = fliplr(AcAll);
    RcAll = fliplr(RcAll);
    figure(1)
    hold on;
    stn = numel(AcAll);
    dx=linspace(-0.4,0.4,stn)*0;
    as = []; 
    rs = [];
    for i = 1:stn % each structure
        scatter(AcAll{i}+dx(i),RcAll{i},[],ones(numel(AcAll{i}),1)*i,'.');
        as = [as; AcAll{i}];
        rs = [rs; RcAll{i}];
        %pause(0.1)
    end


    [N,edges]=histcounts(as,nhc);
    rmn = [];
    asv = [];
    for i = 1:numel(edges)-1
        rmn(i) = mean(rs(  find((edges(i)<as) .* (as<edges(i+1))) ));
        asv(i) = (edges(i)+edges(i+1) )/2;
    end
    plot(asv,rmn,'k')

    hold off;
    colormap('jet')


    set(gcf,'color','w');
    xlabel('intensity')
    ylabel('recruitment number')
    maximize
    ylim([0 0.2]);

    imgFig = getframe(gcf); 
    imgOut = imgFig.cdata;
    imwrite(imgOut,pxStatScatALLFN,'Compression', 'none') 

    save(pxStatALLFN,'RcAll','AcAll')
end

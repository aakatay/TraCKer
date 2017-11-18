%% finds edge and distance to the edge with reference and plot
clear
close all
dbg =1;
fnameRecDIR = rdir('..\*binImgRcrtSum_time*');
fnameADIR = rdir('..\*AVG_lap*'); % pre bleach image
fnameMapDIR=rdir('structMapUpd*.mat');


%isGenStructRec = 1; % generates structures to compare
%%
for j = 1:numel(fnameMapDIR)
    close all;
    fnameMap = fnameMapDIR(j).name;
    cellLabel = fnameMap(13:end-4);
    fnameRec = fnameRecDIR(j).name;
    fnameA = fnameADIR(j).name;
    fnameD = [ 'edgeDistHist' fnameMap(13:end-3) 'tif']; % distribution
    fnameMAT = [ 'edgeDistData' fnameMap(13:end)]; % distribution
    delete(fnameD);
    load(fnameMap); % Lt&Bt
    R = double(imread(fnameRec));
    
    % initialize output arrays
    mxEdgeDist = 10;
    edgeDistRec = nan(1,mxEdgeDist); % recruitment data
    edgeDistRef = nan(1,mxEdgeDist); % reference (uniform)
    edgeDistStr = nan(1,mxEdgeDist); % pre bleach structure
    
    % process all
    Lt = LtNew;
    Bt = BtNew;
    clear LtNew BtNew;
    if 0
        % delete non body structs
        ixdel = find(ps(:,2)~=11);
        if numel(ixdel) == size(LtNew,3) % selection not done yet
            if strcmp(input('select all processed structs(y)?','s'),'y')
                ps(ps(:,2)==1,2)=11;
                save(fnameMap,'BtNew','LtNew','ps','tightPos','-v7.3'); % structMapa
                ixdel = find(ps(:,2)~=11);
            end
        end

        ixdel=[]; 
        ps(ixdel,:)=[];
        Lt(:,:,ixdel) = [];
        Bt(ixdel) = [];
    end
    
    
    if 0 % debug
        ixd=13;
        ps = ps(ixd,:);
        Lt = Lt(:,:,ixd);
        Bt = Bt(ixd);
    end
    
    
    ilast = 0;
    ilastr = 0;
    ilastrA = 0;
    d=[];dr=[];drA=[];dcum=[];drcum=[]; drAcum=[];
        
    hw = waitbar(0,'edge dist calc...');
    setfig=0;
    for k = 1:length(Bt) % each tight section
        %if k~=5,continue; end;
        boundaryT = Bt{k};
        %boundary = boundaryT;
        isdrAmin=1;
        
        szXY = size(R);
        [boundaryT]=extendBoundaries(boundaryT,szXY);
        
        genStructRec; %generates structures 'dr' to compare

        [Y,X,N]=findPos(R.*Lt(:,:,k)); 
        if dbg
            imagesc(R.*Lt(:,:,k))
            hold on;
            plot(boundaryT(:,2),boundaryT(:,1),'r');
            hold off;
        end 
        
        dmin=[];
        j1=1;
        for i = 1:numel(Y)
           j2 = j1 + N(i)-1;
           d_ = sqrt((boundaryT(:,1)-Y(i)).^2+(boundaryT(:,2)-X(i)).^2);
           dmin(j1:j2) = min(d_);
           j1 = j2+1;
        end
        bw = 1;
        d = dmin;
        dr = drmin;
        drA = drAmin;
        waitbar(k/length(Bt),hw,'edge dist calc...');

        %% plot
        figure(2); % distance to the edge plot
        % 1/3 data
        bw=1; % binwidth
        [Nd,xed]=histcounts(d,'BinWidth',bw); % xe: edges
        xc = xed(1:end-1)+bw/2; % center pos
        plot(xc-bw/2,Nd);
        edgeDistRec(k,1:numel(Nd)) = Nd;
        Ns = sum(Nd);
        xm = min(xed);

        hold on
        % 2/3 ref flat
        [Nr,xe]=histcounts(dr,'BinWidth',bw); % Nr :reference
        ixm = find(xe==xm); 
        Nsr = sum(Nr(ixm:end));
        Nr = Ns/Nsr*Nr;
        edgeDistRef(k,1:numel(Nr))=Nr;
        xcr = xe(1:end-1)+bw/2; % center pos
        plot(xcr-bw/2,Nr,'r')

        % 3/3 ref preBleach image
        [NrA,xe]=histcounts(drA,'BinWidth',bw); % Nr :reference
        ixm = find(xe==xm); 
        NsrA = sum(NrA(ixm:end));
        NrA = Ns/NsrA*NrA;
        edgeDistStr(k,1:numel(NrA)) = NrA;
        xcrA = xe(1:end-1)+bw/2; % center pos
        plot(xcrA-bw/2,NrA,'k')
        hold off

        xlabel('distance to the edge')
        ylabel('#pxs')
        title('distribution of distance to the edge')
        grid
        legend('data','ref','image')
        set(gcf,'units','pixels','Position',[800,200,510,550]); 
        set(gca,'units','pixels','Position',[50,50,450,450]);
        imgFig = getframe(gcf);
        dataImg = imgFig.cdata; 
        imwrite(uint16(dataImg),fnameD,'tiff','WriteMode','append','Compression', 'none');        
        dcum = [dcum d];
        drcum = [drcum dr];
        drAcum = [drAcum drA];
        edgeDistMeanArr(k,:) = [sum(d)/numel(d) sum(dr)/numel(dr) sum(drA)/numel(drA) stSz];
    end
    d = dcum;
    dr = drcum;
    drA = drAcum;
    close(hw)
    edgeDistMeanCum = [sum(d)/numel(d) sum(dr)/numel(dr) sum(drA)/numel(drA)]; % [rec ref str]
    edgeDistSumCum = [sum(d) sum(dr)*numel(d)/numel(dr) sum(drA)*numel(d)/numel(drA)]; % [rec ref str]
    save(fnameMAT,'edgeDistRec','edgeDistRef','edgeDistStr','edgeDistMeanCum','edgeDistSumCum','edgeDistMeanArr');
    disp(fnameMAT(13:end-4));
    % edgeDistRatios: print to txt
    fid=fopen('edgeDistRatios.txt','wt');
    fprintf(fid,'mean: \n');
    fprintf(fid,'%.04f\t',edgeDistMeanCum);
    fprintf(fid,'\nsum: \n');
    fprintf(fid,'%.04f\t',edgeDistSumCum);
    fclose(fid);
    
    %% plot scatter maps : dist vs area
    fnameMapDIR=rdir('structMapUpd*.mat');
    cellLabel = fnameMapDIR.name(13:end-4);
    matFNdir = rdir('**\edgeDistData*.mat');
    load(matFNdir.name,'edgeDistMeanArr');
    x = edgeDistMeanArr(:,4); % area
    y = edgeDistMeanArr(:,1)./edgeDistMeanArr(:,2); % area
    scatter(x,y,'.')
    mxx = max(x);
    mxy = max(y);
    nAbove = sum(y>1);
    nBelow = sum(y<=1);
    N = sum(y>0);
    ratAbove = nAbove/N;
    ratBelow = nBelow/N;
    line([0 mxx], [1 1]);
    text(mxx*0.8,mxy*0.9,sprintf('>1:\t%.02f (%i)\n<=1:\t%.02f (%i)',ratAbove,nAbove,ratBelow,nBelow));
    xlabel('area');
    ylabel('edge distance (norm)')
    grid minor
    title(cellLabel)
    
    fnameEdgeDist=sprintf('edgeDistvsArea%s.tif',cellLabel);
    imgFig = getframe(gcf);
    dataImg = imgFig.cdata; 
    imwrite(uint16(dataImg),fnameEdgeDist); % structMap
    
    
    
    
end

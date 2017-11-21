% reads the 'filename_001.tif' and data from tracking of 'filename_002.tif' 
% processes 'filename_001.tif'
function rtTrackSNR

cd('E:\MATLAB\TIRFcalibration\data\Ata01_5_125X100Y50x50_realtime')    
F = findall(0,'type','figure'); delete(F);
    
    cd waSeq\tracker\rtData
    dbgSel = 0;
    dbgSNRimg = 0; % SNR movie
    dbgSNR = 0; % intensity traces
    SNRmul = 1000; % for intesity scale 'SNRimg.tif'
    nt0 = 3000; % # of traces to initialize
    
    
    %% load cfgRT
    cfg = '..\..\..\cfgRT';
    c_=load(cfg);
    cfg = c_.cfg;
    
    ndigit = cfg.ndigit; % # of digits for sequence number
    label = cfg.label;
    w = cfg.w;
    h = cfg.h;
    wsz = cfg.wsz; % window size for SM crop
    szXY = [w h];
    szYX = fliplr(szXY);
    acqTime = cfg.acqTime; % [s]
    stdWin = cfg.stdWin; % number of frames to calc. std

    digitFormat = sprintf('%%0%1ii',ndigit);
    s = wsz/2; 
    mag = 4; % display size
    [mag, pos, szx, szy ] = calcMaxMag(zeros(szXY),mag);
    szXYmag = [szx, szy];
    szYXmag = fliplr(szXYmag);
    SNRmovVoronoi = zeros([szYXmag 3]);
    
    %% input files
    traceFN = 'traceData_';
    fn_ = load('..\..\..\fname0'); % filename_WA_
    fname = fn_.fname0;
    fname = ['..\..\..\' fname];   

    %% output files
    SNRmovieFN = 'SNRmovie.tif';
    SNRmovieConvFN = 'SNRmovieConv.tif';
    SNRmovieVoronoiFN = 'SNRmovieVoronoi.tif';
    SNRplotFN = 'SNRplot.tif';
    numSMplotFN = 'SNRnumSMplot.tif';
    snrSTDplotFN = 'SNRstdplot.tif';
    SMdataFN = 'SMdata.mat';
    
    %% SNR image figure
    if dbgSNRimg
        tit = 'SNR image';
        pos(1) = pos(1) - 800;
        pos(1) = pos(1)- 700;            
        figSNRimg = figure('DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),'Position',[pos/2 szXYmag(1) szXYmag(2)]);
        axeSNR = axes('Parent',figSNRimg,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','Nextplot','replacechildren','XLim',0.5+[0 szXYmag(1)],'YLim',0.5+[0 szXYmag(2)]);
    end

    %% crop
    xc1 = 160;
    yc1 = szXY(2)-110;
    xc2 = xc1+20;
    yc2 = yc1+25;

    xc1 = 1;
    yc1 = szXY(2)-1;
    xc2 = xc1+100;
    yc2 = yc1+100;

    xc1 = 1;
    yc1 = 1;
    xc2 = xc1+50;
    yc2 = yc1+50;
    
    %% feedback to calling function
    fcall = 'rtTrackSNR';
    fdbck.inWait = 0;
    fdbck.inWaitCounting = 0;
    fdbck.inPause = 0;
    fdbck.inSave = 0;
    fdbck.inSaveCounting = 0;
    fdbck.inSaveCountingIX = 0;
    fdbck.inSaveCountingMAX = cfg.inSaveCountingMAX;
    fdbck.inStop = 0;    
    
    
    %% display figure
    if dbgSNR
        figSNR = figure(1212);
        dispSz = [1400 900];
        tit = 'Trace Plots';
        set(figSNR,'DoubleBuffer','on','Menubar','none','Name',tit,'NumberTitle','off','Colormap',gray(256),...
            'Position',[250 50 dispSz(1) dispSz(2)]);
    end
    
    %% first frame
    fnameSeq = [fname num2str(1,digitFormat) '.tif'];
    Afrst = imread(fnameSeq);
    Afrst = padarray(Afrst,[s s]);
    
    %% loop
    n = 2;
    f = 1; % a index
    aIX = cell(nt0,1);
    fIX = cell(nt0,1);
    m = 1; % SM index for SNR
    ixFrm = 1;
    Xs = []; Ys = []; Frm = []; TRinf = []; B = []; INT = []; SNR = []; MIX = [];
    while (1)        
        %% load data
        while (1) % wait for update
            traceSeq = [traceFN label '_' num2str(n,digitFormat) '.mat'];
            if ~exist(traceSeq)
                runAnalysisWin
                return
                [fdbck] = funcFeedback(cfg.msgTXT,fdbck,fcall);
                if fdbck.inStop, break;  end % STOP
            else
                break; % continue
            end 
        end      
        
        fn_ = load(traceSeq); % trace MAT
        TraceX = fn_.TraceX;
        TraceY = fn_.TraceY;
        trInf = fn_.trInf;

        %% read trace data
        X = TraceX;
        Y = TraceY;
        
        % SM selection
        %trInf = trInf(515:518,:);
        %frm1 = trInf(:,1);
        %trInf(frm1~=1)=[]; % remove new recruitments

        frm1 = trInf(:,1);
        frm2 = trInf(:,2)+trInf(:,1)-1;
        nt = size(trInf,1); % # of single molecule
        if size(aIX)<nt, aIX{nt} = []; end % increase array size
        if size(fIX)<nt, fIX{nt} = []; end

        dbg = 0;
        if dbg    
            dbgImg = zeros(szXYmag);
            xdbg = round(trInf(:,4));
            ydbg = round(trInf(:,5));
            dbgImg(sub2ind(size(dbgImg),ydbg,xdbg)) = 1;
            figure(209); imagesc(dbgImg); axis image
            imwrite(dbgImg,'dbgImg.tif')
        end        

        %% PROCESS FRAME ===========================================================
        if dbgSNRimg, snrIMG = zeros(szXYmag); end
        
        fnameSeq = [fname num2str(n,digitFormat) '.tif'];
        A = imread(fnameSeq);
        A = padarray(A,[s s]);
        IX = find((frm1<=n) .* (n<=frm2));
        fTrace = n - frm1(IX) + 1;
        
        x = ceil(cellfun(@(v) v(end), X(IX)));
        y = ceil(cellfun(@(v) v(end), Y(IX)));    
        
        %% select isolated SM
        findSMisolated; % (x,y) -> ixs
        ixSM = IX(ixs);
        if isempty(ixSM), continue; end

        %% crop single molecule windows
        
        tt =0;
        % construct crop image stack
        for j = 1:numel(ixSM) % each single molecule in that frame
            ixsm = ixSM(j); % single molecule index
            y_ = y(j);
            x_ = x(j);
            if ((x_<xc1) || (y_<yc1) || (x_>xc2) || (y_>yc2))
                tt = tt + 1;
                %continue; 
            end
            if ((x_<=s) || (y_<=s) || (x_+s-1>szXY(1)) || (y_+s-1>szXY(2))), continue; end
            if n == 2 % add the first frame data
                a(:,:,f) = Afrst(y_-s:y_+s-1,x_-s:x_+s-1); % crop
                aIX{ixsm} = [aIX{ixsm} f]; % index for the frame of the single molecule
                fIX{ixsm} = [fIX{ixsm} 1]; % frames
                f=f+1;
            end
            a(:,:,f) = A(y_-s:y_+s-1,x_-s:x_+s-1); % crop
            aIX{ixsm} = [aIX{ixsm} f]; % index for the frame of the single molecule
            fIX{ixsm} = [fIX{ixsm} n]; % frames
            f=f+1;
        end
        
        %% calculate SNR
        if n >= stdWin 
            aIXix_ = find(~cellfun(@isempty,aIX)); % all SM found so far
            aIXix = aIXix_;
            aIXix( (   (n-frm1(aIXix)+1)<stdWin   )) = []; % remove SM shorter than stdWin
            aIXix( (   frm2(aIXix) < n  )) = []; % remove SM bleach before the current frame (n)
            mlast = m;
            if dbgSNRimg, snrImg = zeros(szXYmag); end
            num_ixImg2(n) = 0;
            for j = 1:numel(aIXix) % for each SM 
                ix = aIXix(j);
                fr = fIX{ix}; % frames
                if fr(end) ~= n, 
                    if dbgSel, disp('1nodata'); end
                    continue; end % no data in this frame
                if numel(fr)-stdWin+1 < 1,if dbgSel,disp('2notenough');end; continue; end % not enough data
                if fr(end)-stdWin+1 ~= fr(end-stdWin+1) % missing data in the last stdWin frames
                    if dbgSel,disp('3missingdata'); end
                    continue;
                end
                num_ixImg2(n) = num_ixImg2(n) + 1;
                aix = aIX{ix};
                asm = a(:,:,aix); % SM crop stack
                smLast = asm(:,:,end); % last image

                intWin = double(smLast(s-1:s+2,s-1:s+2)); % 4by4
                int = sum(intWin(:));
                peak = max(intWin(:));
                asm(s-1:s+2,s-1:s+2,:) = nan;
                a2D = double(reshape(asm,wsz^2,size(asm,3)));

                astd = std(a2D,0,1);
                astd(isnan(astd)) = [];
                b = mean(astd); % background std
                snr = peak/sqrt(peak+b^2);

                Xs(m) = ceil(cellfun(@(v) v(end), X(ix)));
                Ys(m) = ceil(cellfun(@(v) v(end), Y(ix)));
                Frm(m) = n;
                TRinf(m,:) = trInf(ix,:);     
                B(m) = b;
                INT(m) = int;
                SNR(m) = snr;                   
                MIX(m) = ix; % SM index
                
                if dbgSNRimg
                    Xmag = round(Xs(m)*mag);
                    Ymag = round(Ys(m)*mag);
                    snrImg(Ymag,Xmag,j) = SNR(m);
                end

                m = m + 1;
                %BCK(n) = mean(bckgrnd(:));
                %INT(n) = max(peak(:));
                %RECT(n,:) = rect;
            end
            
            if dbgSNRimg, snrIMG = sum(snrImg,3); end
            %sn{n} = find(snrImg~=0);
            %% display
            if dbgSNR
                figure(figSNR)
                runIntensityTracesRT
            end
            %% save
            save(SMdataFN,'Xs','Ys','Frm','TRinf','B','INT','SNR','MIX','aIX','fIX');
        end
        
        if dbgSNRimg
            figure(figSNRimg);
            set(axeSNR,'Parent',figSNRimg,'DataAspectRatio',[1 1 1],'Position',[0 0 1 1],'Visible','off','XLim',1+[0 szXYmag(1)],'YLim',1+[0 szXYmag(2)]);
            imagesc(flipud(snrIMG))
            pause(acqTime)
            SNRimg(:,:,n) = snrIMG;
        end
        
        %sum(reshape(SNRimg,260^2,226),1)

        if ~isempty(Xs)
            XYS = [Xs' Ys' SNR']; 
            
        end
        ixFrm = [ixFrm m]; % # of spots in each frame

        n = n + 1;
        continue
         
%% STOP

        if dbg
            ix = find(~cellfun(@isempty,aIX));
            ix(1:numel(MIX),2)=MIX';

        end

        if 0  % write SM stack
            %%
            ix = find(~cellfun(@isempty,aIX));
            ixs=ix(1);
            frm = aIX{ixs};
            a(:,:,frm);
            xs = round(trInf(ixs,4));
            ys = round(trInf(ixs,5));
            fname = sprintf('SM%03i_x%iy%i',ixs,xs,ys);
            stackWrite(a,fname)

        end
    end
    
    
    
%% FUNCTIONS --------------------------------------------------------    
    function findSMisolated 
        % find non-overlapping SM
        smMap = zeros(szYX);
        smMap(sub2ind(szYX,y,x)) = 1;
        smMapCv = conv2(smMap,ones(wsz),'same');
        smMapCv(smMapCv~=1) = -10000;
        smMapCv2 = conv2(smMapCv,ones(wsz),'same');
        smMapCv3 = zeros(szYX);
        smMapCv3(smMapCv2 == wsz^2) = 1;
        smMapCv4 = circshift(smMapCv3,[1 1]);
        smMapCv5 = smMapCv4+smMap;
        ixImg = find(smMapCv5==2);
        ixIMG = sub2ind(szYX,y,x)';
        ixs = find(ismember(ixIMG,ixImg));
        

        dbg2 = 0;
        if dbg2 
            figure(344);imagesc(smMap)
            figure(345);imagesc(smMapCv)
            figure(346);imagesc(smMapCv2)
            figure(347);imagesc(smMapCv3)
            figure(348);imagesc(smMapCv4)
            figure(349);imagesc(smMapCv5)
            cc=4;
        end
    end
        
    function runSNRmov
    end

    function runSNRplot
        %% SNR plot
        figSNRplot = figure(201);
        SNR2d = reshape(SNRimg,[szXYmag(2)*szXYmag(1) nfr]);
        SNR2dSUM = sum(SNR2d,1);
        SNR2dNUM = sum(SNR2d>0,1);
        snrMean = SNR2dSUM./SNR2dNUM;

        plot(snrMean); 
        title('SNR mean')
        xlabel('frames')
        ylabel('SNR')   
        imgFig = getframe(gcf);
        imgOut = imgFig.cdata;
        imwrite(imgOut,SNRplotFN);
    end
    
    function runSNRstdPlot
        %% SNR std plot
        figSNRstdPlot = figure(202);
        SNR2d = reshape(SNRimg,[szXYmag(2)*szXYmag(1) nfr]);
        SNR2dSTD = std(SNR2d,1);

        plot(SNR2dSTD); 
        title('SNR STD')
        xlabel('frames')
        ylabel('SNR STD')   
        imgFig = getframe(gcf);
        imgOut = imgFig.cdata;
        imwrite(imgOut,snrSTDplotFN);
    end

    function runNumSMplot
        %% number of SM plot
        figNumSMplot = figure(204);
        plot(num_ixImg2);
        title('number of SM for SNR reading')
        xlabel('frames')
        ylabel('Number of Single Molecules')   
        imgFig = getframe(gcf);
        imgOut = imgFig.cdata;
        imwrite(imgOut,numSMplotFN);
    end    

    function runIntensityTracesRT
        %% intensity traces
        hold on;
        %load(SMdataFN); % 'Xs','Ys','Frm','TRinf','B','INT','SNR','MIX');
        smix = unique(MIX); % # selected SM
        if isempty(smix), return; end
        if smix(1)==0, smix(1)=[]; end % remove zero
        ns = numel(smix);

        intShft = 10e4;
        intShft = 0;
        ixs = [];
        for i = 1:ns
            ix = smix(i); % SM index
            mix = find(MIX==ix);
            frm = Frm(mix);
            FRM(i)=numel(frm);
            %if numel(frm)<217,continue; end
            ixs = [ixs ix];
            int = INT(mix);
            INTsave(:,i) = int;
            plot(frm,int+intShft*(i-1));
        end
        save('SNRintTraces','INTsave','frm','ixs')
        savefig('SNRintTraces');
        hold off;
    end

    function runIntensityTraces
        %% intensity traces
        figIntensityTraces = figure(1212); maximize;clf;
        hold on;
        %load(SMdataFN); % 'Xs','Ys','Frm','TRinf','B','INT','SNR','MIX');
        smix = unique(MIX); % # selected SM
        if smix(1)==0, smix(1)=[]; end % remove zero

        smixSel = [11 14 16 21];
        ns = numel(smix);
        ns2 = 0;
        ns3 = 0;
        FRM = [];
        intShft = 10e4;
        intShft = 0;
        ixs = [];
        for i = 1:ns
            ixc = smix(i); % SM index
            mix = find(MIX==ixc);
            frm = Frm(mix);
            FRM(i)=numel(frm);
            if numel(frm)<100,continue; end
            if numel(frm)<217,continue; end

            ns2 = ns2 + 1;
            %if ~ismember(ns2,smixSel),continue; end;
            ns3= ns3+1;
            ixs = [ixs ixc];
            int = INT(mix);
            INTsave(:,ns3) = int;
            subplot(4,1,ns3)
            plot(frm,int+intShft*(ns3-1));
            ylim([0 max(int)*1.2])
            pause(0.5)
        end
        save('SNRintTraces','INTsave','frm','ixs')
        savefig('SNRintTraces');
        hold off;
    end


%% ============ Analysis ===================================================
    function runAnalysisWin
        fig= uifigure('Position',[20 200 300 300]);

        x1 = 20;
        dy = 30;
        y1 = [0:3]*dy+10;

        w = 100;
        h1 = 20;
        btn1 = uibutton(fig,'Position',[x1 y1(1) w h1],'Text','Intensity','ButtonPushedFcn', @(btn1,event) dispTraceInt);
        btn2 = uibutton(fig,'Position',[x1 y1(2) w h1],'Text','SNR','ButtonPushedFcn', @(btn2,event) runPause);
        btn3 = uibutton(fig,'Position',[x1 y1(3) w h1],'Text','SNR movie','ButtonPushedFcn', @(btn3,event) writeSNRmov);
        btn4 = uibutton(fig,'Position',[x1 y1(4) w h1],'Text','SNR voronoi mov','ButtonPushedFcn', @(btn4,event) writeSNRvoronoiMov);
        while(1)
            pause(1)
        end
    end

    function dispTraceInt
        figure(901)
        hold on
        smix = unique(MIX); % # selected SM
        if isempty(smix), return; end
        if smix(1)==0, smix(1)=[]; end % remove zero
        ns = numel(smix);

        intShft = 10e4;
        intShft = 0;
        ixs = [];
        for i = 1:ns % each selected trace
            ix = smix(i); % SM index
            mix = find(MIX==ix);
            frm = Frm(mix);
            FRM(i)=numel(frm);
            %if numel(frm)<217,continue; end
            ixs = [ixs ix];
            int = INT(mix);
            %INTsave(:,i) = int;
            plot(frm,int+intShft*(i-1));
        end
        save('SNRintTraces','INTsave','frm','ixs')
        savefig('SNRintTraces');      
        hold off  
        
    end


    function writeSNRvoronoiMov
        delete(SNRmovieVoronoiFN);
        
        CM = jet(256);
        CM = gray(256);
        for i = 1:n-2
            xys = XYS(ixFrm(i):ixFrm(i+1)-1,:);
            if size(xys,1)<3, continue; end
            snrDisp = xys(:,3);
            if max(snrDisp)<1.5
                sscale = 1.5;
                EdgeColorSel = 1; % red
            elseif max(snrDisp)<5
                sscale = 5;
                EdgeColorSel = 2; % green 
            else % max(snrDisp)<10
                sscale = 10;
                EdgeColorSel = 3; % blue
            end
            SNRmovVoronoi = getVoronoinImg(xys(:,1:2),snrDisp/sscale,szXY,mag,CM,EdgeColorSel);
            imwrite(SNRmovVoronoi,SNRmovieVoronoiFN,'WriteMode','append','Compression', 'none') 
        end
    end

    function writeSNRmov
        if ~dbgSNRimg, return; end
        %% SNR movie
        if ~isempty(find(SNRimg*SNRmul>=2^16)), warning('SNR image saturated use lower a SNRmul');end
        SNRmov = SNRimg*SNRmul;
        stackWrite(SNRmov,SNRmovieFN);
        cvWin = ones(5);
        SNRmovCv = convn(SNRmov,cvWin,'same');
        stackWrite(SNRmovCv,SNRmovieConvFN);
    end

end































function rtTraCKerTrace
% minTraceLength = 2
% traceJmpForCombination = 1
%cd('E:\MATLAB\TIRFcalibration\data\Ata01_5_125X100Y50x50_realtime'); 
    %cd waSeq\tracker\
    %% files
    
    cfg_ = load('..\..\cfgRT');
    cfg = cfg_.cfg;
    
    
    
    
    
    
    sptJmpForTracing = cfg.sptJmpForTracing;
    dn = cfg.sptReAppearTime;
    ndigit = cfg.ndigit; % # of digits for sequence number
    label = cfg.label; % # of digits for sequence number
    w = cfg.w;
    h = cfg.h;
    szXY = [w h];
    szYX = fliplr(szXY);
        wsz = cfg.wszTracker; % window size for SM crop
    %wsz = 5; % window size for SM crop
    w = floor(wsz/2);
    outDIR = 'rtData\';
    
    tic; logFN = cfg.logTrace; fid = fopen(logFN,'w'); wait = 0;
    clck = clock; fprintf(fid,'start time m= %2i secs=%6.03f\n',clck(5),clck(6));
    
    % prime numbers
    if dn>2, error('update the code (sptReAppearTime)'); end
    pNlim = [200000 400000 600000 800000];
    pns = genPrimeNumSet(pNlim); % pNlim --> pns : (3,10000)
            
    NT1_ = 0;
    %%
    digitFormat = sprintf('''%%0%1ii''',ndigit);
    

    %% feedback to calling function
    fcall = 'rtTraCKerTrace';
    fdbck.inWait = 0;
    fdbck.inWaitCounting = 0;
    fdbck.inPause = 0;
    fdbck.inSave = 0;
    fdbck.inSaveCounting = 0;
    fdbck.inSaveCountingIX = 0;
    fdbck.inSaveCountingMAX = cfg.inSaveCountingMAX;
    fdbck.inStop = 0;    
    
    %% loop
    ixSptFrm = 1; % first spot in the frame
    n = 1; % frame number
    ixTr = 1; % trace number
    X=[];Y=[];INT=[];
    TraceX={};TraceY={};
    trInf = [];
    pnIMG = nan([szYX dn+1]);
    while (1)
        time = toc; fprintf(fid,'while loop n=%3i time=%6.03f\n',n,time);
        %% output filename
        digitTXT = eval(['sprintf(' digitFormat ',n)'] );
        traceDataFileNm = [outDIR 'traceData_' label '_' digitTXT '.mat'];
        
        %% load data
        while (1) % wait for update
            tloop(n)=toc;
            posFN = rdir(['posData-coeff*_' label '_' digitTXT '.mat']);
            if isempty(posFN)
                if wait == 0, time = toc; fprintf(fid,'wait for   n=%3i time=%6.03f\n',n,time); wait = 1; end
                [fdbck] = funcFeedback(cfg.msgTXT,fdbck,fcall);
                if fdbck.inStop, break;  end % STOP
            else
                time = toc; wait = 0; fprintf(fid,'updated    n=%3i time=%6.03f\n',n,time);
                break; % continue
            end 
            if 1 
                TLOOP = tloop(2:end)-tloop(1:end-1);
                nv = 1:n-1;
                %plotyy(nv,TLOOP,nv,tsave)
                plot([TLOOP;tsave;tprm]')
                title(sprintf('mean time : %.03f',mean(tprm)));
                legend('loop' ,'save','primes')
                return;
            end
            %pause(0.010)
        end
        posfn = posFN.name;
        
        posfn_=load(posfn);
        XC = posfn_.X;
        YC = posfn_.Y;
        INTC = posfn_.INT;
        INT = [INT INTC];
        npos = numel(XC); % # of localizations
        
        snc = rem(n,3); % prime number set number (current)
        snp = rem(n-1,3); % prime number set number (previous)
        sng = rem(n-2,3); % prime number set number (gap: 2frm before)
        if snc==0,snc=3;elseif snp==0,snp=3;else,sng=3;end
        
        pnImg = genPrimeNumImg(YC,XC,szYX,pns,snc,wsz,npos);
        pnIMG(:,:,1) = pnImg;
        pnIMG1 = pnIMG(:,:,1);
        %pnIMG1 = abs(pnIMG(:,:,2));
        uniqImg = unique(pnIMG1(~isnan(pnIMG1(:))));
        uprm = isprime(uniqImg);
        uprm(uprm>pns(snc,numel(YC)))=0;
        
        nt1 = sum(uprm);
        ixtc = nan(1,nt1); % trace indices in the current frame
        
        if 0         
            pnIMG2 = abs(pnIMG(:,:,2));
            prms = unique(pnIMG2(~isnan(pnIMG2)));
            mxPrm = max(prms(isprime(prms)));
            nt2 = sum(isprime(prms));
        end
        
        nSPOTs(n) = length(XC); % number of spots in that frame
        %ixSpt = ixSptFrm(n):ixSptFrm(n)+nSPOTs(n)-1; % indices of the spots in the current frame
        %frmNoSpot(ixSpt) = n;
        ixSptFrm(n+1) = ixSptFrm(n)+nSPOTs(n); % first spot in the frame
        if n == 1 % first frame only
            n = 2; 
            XP = XC; YP = YC;
            ixtp = ixtc;
            pnIMG = circshift(pnIMG,[0 0 1]);
            continue; 
        end
        
        %% 1: check prev frame
        pnMATCH = pnIMG(:,:,1).*pnIMG(:,:,2); 
        %%
        pnMatch = unique(pnMATCH(~isnan(pnMATCH)));
        nt = numel(unique(abs(pnMatch)));
        if numel(pnMatch) ~= nt % BUG
            ccc=2;
        end
        if 0 
            pnIMG2 = pnIMG(:,:,2);
            nt2 = sum(isprime(unique(pnIMG2(~isnan(pnIMG2)))));
            figure(11); imagesc(pnIMG(:,:,1)); figure(12); imagesc(pnIMG(:,:,2)); 
        end
        tprm(n) = 0;
        for i = 1:nt % each trace
            if pnMatch(i) == 0, continue;end % no matching with prev frame
            %fct = factor(abs(pnMatch(i)));
            tprm0 = toc;
            fct = factorFast(abs(pnMatch(i)),pns(snc,:),pns(snp,:));
            tprm2 = toc;
            tprm(n) = tprm(n)+tprm2-tprm0;
            if numel(fct) ~=2, a='disp(fct)';continue; % overlapping traces
            else % match found : (numel(fct)=2)
                [fixc,fixp] = getFix(n,'c','p');
                ixc = find(fct(fixc) == pns(snc,:)); % index to XC YC arrays
                ixp = find(fct(fixp) == pns(snp,:)); % index to XP YP arrays
                pnIMG(:,:,1) = pnIMGremoveLocalization(pnIMG(:,:,1),YC(ixc),XC(ixc),w,szXY);
                
                %tix_ = find(ismember(ixtp,ixp)); % check traces in prev frame
                tix = ixtp(ixp); % check traces in prev frame
                if  ~isnan(tix) % add to the trace
                    %tix = ixtp(tix_);
                    
                    TraceX{tix} = [TraceX{tix} XC(ixc)];
                    checkTX(1)
                    TraceY{tix} = [TraceY{tix} YC(ixc)];
                    trInf(tix,2) = trInf(tix,2)+1; % num frames
                    
                    ixtc(ixc) = tix; % index to TraceX 
                    ixtp(ixtp==tix)=nan;
                elseif pnMatch(i)>0 % new trace
                    tracex = [XP(ixp) XC(ixc)];
                    tracey = [YP(ixp) YC(ixc)];
                    TraceX{end+1} = tracex;
                    checkTX(2)
                    TraceY{end+1} = tracey;
                    trInf(end+1,1) = n-1; % frst frame
                    trInf(end,2) = 2; % num frames
                    ixtc(ixc) = numel(TraceX); % index to TraceX 
                    pnIMG(:,:,2) = pnIMGremoveLocalization(pnIMG(:,:,2),YP(ixp),XP(ixp),w,szXY);
                end
            end
        end
        
        %% 2: check 2frm before (with a gap)
        if n == 2 % second frame only
            n = 3; 
            XG = XP; YG = YP;
            XP = XC; YP = YC;
            ixtg = ixtp;
            ixtp = ixtc;
            pnIMG = circshift(pnIMG,[0 0 1]);
            tsave0 = toc;
            save(traceDataFileNm,'TraceX','TraceY','trInf')
            tsave2 = toc;
            tsave(n) = tsave2-tsave0;
            continue; 
        end
        pnIMGr = pnIMG(:,:,1); % remaining
        pnIMGr(pnIMGr<0)=0;
        pnMATCH = pnIMGr.*pnIMG(:,:,3); 
        pnMatch = unique(pnMATCH(~isnan(pnMATCH)));
        nt = numel(unique(abs(pnMatch)));
        for i = 1:nt % each trace
            if pnMatch(i) == 0, continue;end % no matching with prev frame
            %fct = factor(abs(pnMatch(i)));
            fct = factorFast(abs(pnMatch(i)),pns(snc,:),pns(sng,:));
            if numel(fct) ~=2, a='disp(fct)';continue; % overlapping traces
            elseif numel(fct) == 1, continue; % no matching with prev frame
            else % match found : (numel(fct)=2)
                [fixc,fixg] = getFix(n,'c','g');
                ixc = find(fct(fixc) == pns(snc,:)); % index to XC YC arrays
                ixg = find(fct(fixg) == pns(sng,:)); % index to XP YP arrays
            
                pnIMG(:,:,1) = pnIMGremoveLocalization(pnIMG(:,:,1),YC(ixc),XC(ixc),w,szXY);
                %tix_ = find(ismember(ixtg,ixg)); % check traces in 2prev frame
                tix = ixtg(ixg); % check traces in prev frame
           
                if  ~isnan(tix) % add to the trace
                    xg = (TraceX{tix}(end)+XC(ixc))/2;
                    yg = (TraceY{tix}(end)+YC(ixc))/2;

                    TraceX{tix} = [TraceX{tix} xg XC(ixc)];
                    checkTX(1)
                    TraceY{tix} = [TraceY{tix} yg YC(ixc)];
                    trInf(tix,2) = trInf(tix,2)+2; % num frames
                    ixtc(ixc) = tix; % index to TraceX 
                    ixtg(ixtg==tix)=nan;
                elseif pnMatch(i)>0 % new trace
                    xg = (XG(ixg)+XC(ixc))/2;
                    yg = (YG(ixg)+YC(ixc))/2;
                    tracex = [XG(ixg) xg XC(ixc)];
                    tracey = [YG(ixg) yg YC(ixc)];
                    TraceX{end+1} = tracex;
                    checkTX(2)
                    TraceY{end+1} = tracey;
                    trInf(end+1,1) = n-2; % frst frame
                    trInf(end,2) = 3; % num frames
                    ixtc(ixc) = numel(TraceX); % index to TraceX 
                    pnIMG(:,:,3) = pnIMGremoveLocalization(pnIMG(:,:,3),YG(ixg),XG(ixg),w,szXY);
                end
            end
        end
        pnIMG = circshift(pnIMG,[0 0 1]);
        XG = XP; YG = YP; 
        XP = XC; YP = YC; 
        ixtg = ixtp; 
        ixtp = ixtc;
        
        tsave0 = toc;
        save(traceDataFileNm,'TraceX','TraceY','trInf')
        tsave2 = toc;
        tsave(n) = tsave2-tsave0;
        n = n + 1;
        cc= 3;
    end
           
    % ====================================================================
    function fct = factorFast(num,pns1,pns2)
        fct=[];
        fctIx1 = find(ismember(num./pns1,pns2));
        if ~isempty(fctIx1)
            fct(1) = pns1(fctIx1);
            fct(2) = num/fct(1);
            fct = sort(fct);
            ccc=4;
        end
        
    end

    function checkTX(nin)
        return;
        NT1 = numel(TraceX{1});
        if NT1_ < NT1
            n
            NT1
            cccc=2;
        end
        NT1_ = NT1;
        if nin==1

            if numel(TraceX{tix})>n
                ccccc=3;
            end
        else
            
            if numel(TraceX{end})>n
                ccccc=3;
            end
        end
    end
    
    
    function pnImg = genPrimeNumImg(YC,XC,szYX,pns,snc,wsz,npos)
        %  -- > pnImg
        %sn: set number
        % max 1e4 localizations in a frame
        smMap = zeros(szYX);
        smMap(sub2ind(szYX,ceil(YC),ceil(XC))) = pns(snc,1:npos);
        pnImg = int32(conv2(int32(smMap),ones(wsz),'same'));
    end
    
    function pns = genPrimeNumSet(pNlim)
        pn_0 = primes(pNlim(1));
        
        % set 1
        npn_0 = numel(pn_0);
        pn_1 = primes(pNlim(2));
        pns(1,:) = pn_1(npn_0+1:npn_0+1e4);
        % set 2
        npn_1 = numel(pn_1);
        pn_2 = primes(pNlim(3));
        pns(2,:) = pn_2(npn_1+1:npn_1+1e4);
        % set 3
        npn_2 = numel(pn_2);
        pn_3 = primes(pNlim(4));
        pns(3,:) = pn_3(npn_2+1:npn_2+1e4);
    end

    function pnIMGr = pnIMGremoveLocalization(pnIMGr,yc,xc,w,szXY)
        x1 = round(xc-w);
        x2 = round(xc+w);
        y1 = round(yc-w);
        y2 = round(yc+w);
        if x1<1, x1=1;end; if y1<1, y1=1; end
        if x2>szXY(1), x2=szXY(1);end; if y2>szXY(2), y2=szXY(2); end
        pnIMGr(y1:y2,x1:x2) = -pnIMGr(y1:y2,x1:x2); % remove the current localization
    end
    
    function [fixc,fixpg] = getFix(n,~,frm2)
    % fixes the dynamic orderin between factor sets
        Fix = [1 2 3]; % index to factors
        Fix = circshift(Fix,-n);
        FixC = Fix(3);
        switch frm2
            case 'p'
                FixPG = Fix(2);
            case 'g'
                FixPG = Fix(1);
        end
        if FixC<FixPG
            fixc = 1;
            fixpg = 2;
        else
            fixc = 2;
            fixpg = 1;
        end
    end
end
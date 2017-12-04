function rtTraCKerTrace(varargin)
% minTraceLength = 2
% traceJmpForCombination = 1
%cd('E:\MATLAB\TIRFcalibration\data\Ata01_5_125X100Y50x50_realtime'); 

% FASTER: try cuda for ismember function in factorFast1()

    isCallOutside = 0;
    if nargin == 1
        isCallOutside = 1;
        cfg = varargin{1};
        PWD = cfg.pwd;
        cd(PWD)
        q=0;
        pause(1);
        while q == 0
            try 
                cd waSeq\tracker\
                q = 1;
            catch
                q=0;
                %save PWD PWD
                %pause(1);
            end
        end
    end
    %PWD = pwd;     save PWD PWD
    
    %% files
    ffT =[];
    cfg_ = load('..\..\cfgRT');
    cfg = cfg_.cfg;
    np = cfg.mxNumLocalization; % max localizations in a frame
        
    sptJmpForTracing = cfg.sptJmpForTracing;
    dn = cfg.sptReAppearTime;
    ndigit = cfg.ndigit; % # of digits for sequence number
    label = cfg.label; % # of digits for sequence number
    w = cfg.w;
    h = cfg.h;
    isTlog = cfg.isTlog;
    tic;
    if isTlog, tic; logFN = cfg.logTrace; fid = fopen(logFN,'w'); wait = 0; end
    if isTlog, clck = clock; fprintf(fid,'start time m= %2i secs=%6.03f\n',clck(5),clck(6)); end
    
    szXY = [w h];
    szYX = fliplr(szXY);
        wsz = cfg.wszTracker; % window size for SM crop
    %wsz = 5; % window size for SM crop
    w = floor(wsz/2);
    outDIR = 'rtData\';
    
    
    % prime numbers
    if dn>2, error('update the code (sptReAppearTime)'); end
    pNlim = [200000 400000 600000 800000];
    [pns,pnsMul] = genPrimeNumSet(pNlim); % pNlim --> pns : (3,10000)
            
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
    pnIMG = zeros([szYX dn+1],'double');
    while (1)

tprm0 = toc; % vvvvvvvvvvvvvvvvvvvvvvvvv
        if isTlog, time = toc; fprintf(fid,'while loop n=%3i time=%6.03f\n',n,time); end
        %% output filename
        digitTXT = eval(['sprintf(' digitFormat ',n)'] );
        traceDataFileNm = [outDIR 'traceData_' label '_' digitTXT '.mat'];
        
        %% load data
        while (1) % wait for update
            tloop(n)=toc;
            posFN = rdir(['posData-coeff*_' label '_' digitTXT '.mat']);
            if isempty(posFN)
                if isTlog, if wait == 0, time = toc; fprintf(fid,'wait for   n=%3i time=%6.03f\n',n,time); wait = 1; end; end
                [fdbck] = funcFeedback(cfg.msgTXT,fdbck,fcall);
                if fdbck.inStop, break;  end % STOP
            else
                if isTlog, time = toc; wait = 0; fprintf(fid,'updated    n=%3i time=%6.03f\n',n,time); end
                break; % continue
            end 
            if exist('tprm')
                save('rtTraCKerTraceTIME.mat','tprm','tloop','ffT'); 
                if ~isCallOutside
                    TLOOP=[];
                    dbgTraceTIME
                    return;
                end
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
tprm2 = toc;
tprm(n,5) = tprm2-tprm0;
        
tprm0 = toc;% vvvvvvvvvvvvvvvvvvvvvvvvv
        pnImg = genPrimeNumImg(YC,XC,szYX,pns,snc,wsz,npos);
        pnIMG(:,:,1) = pnImg;
        pnIMG1 = pnIMG(:,:,1);
        %pnIMG1 = abs(pnIMG(:,:,2));
        uniqImg = unique(pnIMG1(~isnan(pnIMG1(:))));
        uprm = isprime(uniqImg);
        uprm(uprm>pns(snc,numel(YC)))=0;
        
        nt1 = sum(uprm);
        ixtc = nan(1,nt1); % trace indices in the current frame
tprm2 = toc;
tprm(n,4) = tprm2-tprm0;
        
        if 0         
            pnIMG2 = abs(pnIMG(:,:,2));
            prms = unique(pnIMG2(~isnan(pnIMG2)));
            mxPrm = max(prms(isprime(prms)));
            nt2 = sum(isprime(prms));
        end
        
        if n == 1 % first frame only
            n = 2; 
            XP = XC; YP = YC;
            ixtp = ixtc;
            pnIMG = circshift(pnIMG,[0 0 1]);
            continue; 
        end

tprm0 = toc; % vvvvvvvvvvvvvvvvvvvvvvvvv        
        %% 1: check prev frame
        pnMATCH = pnIMG(:,:,1).*pnIMG(:,:,2); 
        pnMatch = unique(pnMATCH(~isnan(pnMATCH)))';
        nm = numel(pnMatch); % # of matchings
        
        pns1 = pns(snc,:)';
        %PNS1 = repmat(pns1,1,nm);
        pns2 = pns(snp,:);
        pnsMulc = pnsMul(snp,:); % current multiple set
        
        %FCT = factorFast1(abs(pnMatch),PNS1,pns2);
        [FCT,mv] = factorFast1(abs(pnMatch),pns1,pns2,pnsMulc); % mv: match vector
        nt = size(FCT,2);
tprm2 = toc;
tprm(n,1) = tprm2-tprm0;
        
        
        if 0 
            pnIMG2 = pnIMG(:,:,2);
            nt2 = sum(isprime(unique(pnIMG2(~isnan(pnIMG2)))));
            figure(11); imagesc(pnIMG(:,:,1)); figure(12); imagesc(pnIMG(:,:,2)); 
        end
        xcrC = []; ycrC = [];
        xcrP = []; ycrP = [];
tprm01 = toc;% vvvvvvvvvvvvvvvvvvvvvvvvv
        for i = 1:nt % each trace
            mix = mv(i);
            %fct = factor(abs(pnMatch(i)));
            %fct = factorFast(abs(pnMatch(i)),pns(snc,:),pns(snp,:));
            [fixc,fixp] = getFix(n,'c','p');
            ixc = find(FCT(fixc,i) == pns(snc,:)); % index to XC YC arrays
            ixp = find(FCT(fixp,i) == pns(snp,:)); % index to XP YP arrays
            
            xcrC = [xcrC XC(ixc)]; % mark in the image
            ycrC = [ycrC YC(ixc)];

            %tix_ = find(ismember(ixtp,ixp)); % check traces in prev frame
            tix = ixtp(ixp); % check traces in prev frame
            if  ~isnan(tix) % add to the trace
                %tix = ixtp(tix_);

                TraceX{tix} = [TraceX{tix} XC(ixc)];
                %checkTX(1)
                TraceY{tix} = [TraceY{tix} YC(ixc)];
                trInf(tix,2) = trInf(tix,2)+1; % num frames
                
    
                ixtc(ixc) = tix; % index to TraceX 
                ixtp(ixtp==tix)=nan;
            elseif pnMatch(mix)>0 % new trace
                tracex = [XP(ixp) XC(ixc)];
                tracey = [YP(ixp) YC(ixc)];
                TraceX{end+1} = tracex;
                %checkTX(2)
                
                TraceY{end+1} = tracey;
                trInf(end+1,1) = n-1; % frst frame
                trInf(end,2) = 2; % num frames
                ixtc(ixc) = numel(TraceX); % index to TraceX 
                xcrP = [xcrP XP(ixp)]; % mark in the image
                ycrP = [ycrP YP(ixp)];
            end
        end
tprm02 = toc;
tprm(n,2) = tprm02-tprm01;

tprm01 = toc;% vvvvvvvvvvvvvvvvvvvvvvvvv
        % mark in the image
        pnIMG(:,:,1) = pnIMGremoveLocalization(pnIMG(:,:,1),xcrC,ycrC,wsz);
        pnIMG(:,:,2) = pnIMGremoveLocalization(pnIMG(:,:,2),xcrP,ycrP,wsz);
tprm02 = toc;
tprm(n,3) = tprm02-tprm01;
        
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
tprm(n,6) = tsave2-tsave0;
            continue; 
        end

        isgap =0;
        if isgap

            pnIMGr = pnIMG(:,:,1); % remaining
            pnIMGr(pnIMGr<0)=0;
            pnMATCH = pnIMGr.*pnIMG(:,:,3); 
            pnMatch = unique(pnMATCH(~isnan(pnMATCH)))';
            nm = numel(pnMatch);
            mv = 1:nm; % match vector

            pns1 = pns(snc,:)';
            PNS1 = repmat(pns1,1,nm);
            pns2 = pns(sng,:);


            FCT = factorFast2(abs(pnMatch),PNS1,pns2);
            nt = size(FCT,2);

            xcrC = []; ycrC = [];
            xcrG = []; ycrG = [];
            for i = 1:nt % each trace
                mix = mv(i);
                if pnMatch(mix) == 0, continue;end % no matching with prev frame
                %fct = factor(abs(pnMatch(i)));
                %fct = factorFast(abs(pnMatch(i)),pns(snc,:),pns(sng,:));
                [fixc,fixg] = getFix(n,'c','g');
                ixc = find(FCT(fixc,i) == pns(snc,:)); % index to XC YC arrays
                ixg = find(FCT(fixg,i) == pns(sng,:)); % index to XP YP arrays

                xcrC = [xcrC XC(ixc)]; % mark in the image
                ycrC = [ycrC YC(ixc)];

                %tix_ = find(ismember(ixtg,ixg)); % check traces in 2prev frame
                tix = ixtg(ixg); % check traces in prev frame
                if  ~isnan(tix) % add to the trace
                    xg = (TraceX{tix}(end)+XC(ixc))/2;
                    yg = (TraceY{tix}(end)+YC(ixc))/2;

                    TraceX{tix} = [TraceX{tix} xg XC(ixc)];
                    %checkTX(1)
                    TraceY{tix} = [TraceY{tix} yg YC(ixc)];
                    trInf(tix,2) = trInf(tix,2)+2; % num frames
                    ixtc(ixc) = tix; % index to TraceX 
                    ixtg(ixtg==tix)=nan;
                elseif pnMatch(mix)>0 % new trace
                    xg = (XG(ixg)+XC(ixc))/2;
                    yg = (YG(ixg)+YC(ixc))/2;
                    tracex = [XG(ixg) xg XC(ixc)];
                    tracey = [YG(ixg) yg YC(ixc)];
                    TraceX{end+1} = tracex;
                    %checkTX(2)
                    TraceY{end+1} = tracey;
                    trInf(end+1,1) = n-2; % frst frame
                    trInf(end,2) = 3; % num frames
                    ixtc(ixc) = numel(TraceX); % index to TraceX 
                    xcrG = [xcrG XG(ixg)]; % mark in the image
                    ycrG = [ycrG YG(ixg)];
                end
            end
            % mark in the image
            pnIMG(:,:,1) = pnIMGremoveLocalization(pnIMG(:,:,1),xcrC,ycrC,wsz);
            pnIMG(:,:,3) = pnIMGremoveLocalization(pnIMG(:,:,3),xcrG,ycrG,wsz);
        end
        
        
        pnIMG = circshift(pnIMG,[0 0 1]);
        XG = XP; YG = YP; 
        XP = XC; YP = YC; 
        ixtg = ixtp; 
        ixtp = ixtc;
        
tsave0 = toc; % vvvvvvvvvvvvvvvvvvvvvvvvv
        save(traceDataFileNm,'TraceX','TraceY','trInf')
tsave2 = toc;
tsave(n) = tsave2-tsave0;
        n = n + 1;
        cc= 3;
    end
           
    % ====================================================================

    function [FCT,nix] = factorFast1(num,pns1,pns2,pnms)
ffT1 = toc;
        % pnms: prime number multiple set
        if num(1)==0,num(1)=[];end
        pnms2D = reshape(pnms,np,np);
ffT12 = toc;
        %[ism,numix] = ismember(pnms2D,num);
        [~,pix,nix] = intersect(pnms2D,num);
        %[~,pix,nix] = fastintersect(pnms2D,num);
ffT13 = toc;
ffT(n,2) = ffT13-ffT12;
        %numix = unique(sparse(numix));
        %if numix(1)==0,numix(1)=[];end
        %[iy_, ix_] = find(ism);
        [iy, ix] = ind2sub(size(pnms2D),pix);
        %[FCTix1, numIx] = find(ismember(num./pnms));
        FCT(1,:) = pns2(iy);
        FCT(2,:) = pns1(ix);
        FCT = sort(FCT);
        %ixr = FCT(1,:).*FCT(2,:) == num(numIx); FCT = FCT(:,ixr);
        %[(FCT(1,:).*FCT(2,:))'-num(1:28)']
ffT2 = toc;
ffT(n,1) = ffT2-ffT1;
    end

    function [pix,nix] = fastintersect(P,N)
        minP = min(P(:));
        shft = minP - 1;
        Ps = P - shft;
        Ns = N - shft;
        maxPs = max(Ps(:));
        Pix = zeros(1, max(maxPs));
        Nix = Pix;
        Pix(P) = 1;
        Nix(N) = 1;
        nix = logical(Pix(N));
        pix = logical(Nix(P));
       
    end

    function FCT = factorFast2(num,pns1,pns2)
    end

    
    function pnImg = genPrimeNumImg(YC,XC,szYX,pns,snc,wsz,npos)
        %  -- > pnImg
        %sn: set number
        % max np localizations in a frame
        smMap = zeros(szYX,'double');
        smMap(sub2ind(szYX,ceil(YC),ceil(XC))) = pns(snc,1:npos);
        pnImg = double(conv2(smMap,ones(wsz,'double'),'same'));
    end
    
    function [pns,pnsMul] = genPrimeNumSet(pNlim)
        pn_0 = primes(pNlim(1));
        
        % set 1
        npn_0 = numel(pn_0);
        pn_1 = primes(pNlim(2));
        pns(1,:) = pn_1(npn_0+1:npn_0+np);
        % set 2
        npn_1 = numel(pn_1);
        pn_2 = primes(pNlim(3));
        pns(2,:) = pn_2(npn_1+1:npn_1+np);
        % set 3
        npn_2 = numel(pn_2);
        pn_3 = primes(pNlim(4));
        pns(3,:) = pn_3(npn_2+1:npn_2+np);
        
        pns12 = pns(1,:)' * pns(2,:);
        pnsMul = pns12(:)'; clear pns12;
        pns23 = pns(2,:)' * pns(3,:);
        pnsMul(2,:) = pns23(:)'; clear pns23;
        pns31 = pns(3,:)' * pns(1,:);
        pnsMul(3,:) = pns31(:)'; clear pns31;
    end

    function pnIMGr = pnIMGremoveLocalization(pnIMGr,yc,xc,wsz)
        pnimgNAN = zeros(size(pnIMGr),'double');
        pnimgNAN(sub2ind(size(pnIMGr),round(yc),round(xc))) = 1;
        wnan = ones(wsz,'double');
        pnimgNAN = double(conv2(pnimgNAN,wnan,'same'));
        pnimgNAN = pnimgNAN + 1;
        pnimgNAN(pnimgNAN>1) = -1;
        pnIMGr = pnIMGr.*pnimgNAN;
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




























    function FCT = factorFastOld(num,pns1,pns2)
        [FCTix1, numIx] = find(ismember(num./pns1,pns2));
        FCT(1,:) = pns1(FCTix1,1);
        FCT(2,:) = num(numIx)./FCT(1,:);
        %ixr = FCT(1,:).*FCT(2,:) == num(numIx); FCT = FCT(:,ixr);
        FCT = sort(FCT);
    end

    function FCT = factorFast1cuda(num,pns1,pns2,pnms)
        % pnms: prime number multiple set
        pnms2D = reshape(pnms,np,np);
        try
            parfor tt=1:numel(num)
                [iyy,ixx] = find(pnms2D == num(tt));
                if ~isempty(iyy)
                    iy(tt) = iyy;
                    ix(tt) = ixx;
                end
            end
        catch ME
            xcxx=2;
        end
        iy(iy==0)=[];
        ix(ix==0)=[];
          
        %[iy, ix] = find(ism);
        %[FCTix1, numIx] = find(ismember(num./pnms));
        FCT(1,:) = pns2(iy);
        FCT(2,:) = pns1(ix);
        FCT = sort(FCT);
        %ixr = FCT(1,:).*FCT(2,:) == num(numIx); FCT = FCT(:,ixr);
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
    

    function [FCT,numix] = factorFast1ismember(num,pns1,pns2,pnms)
ffT1 = toc;
        % pnms: prime number multiple set
        if num(1)==0,num(1)=[];end
        pnms2D = reshape(pnms,np,np);
ffT12 = toc;
        [ism,numix] = ismember(pnms2D,num);
        [numsel,pix,nix] = intersect(pnms2D,num);
ffT13 = toc;
ffT(n,2) = ffT13-ffT12;
        numix = unique(sparse(numix));
        if numix(1)==0,numix(1)=[];end
        [iy, ix] = find(ism);
        %[FCTix1, numIx] = find(ismember(num./pnms));
        FCT(1,:) = pns2(iy);
        FCT(2,:) = pns1(ix);
        FCT = sort(FCT);
        %ixr = FCT(1,:).*FCT(2,:) == num(numIx); FCT = FCT(:,ixr);
        %[(FCT(1,:).*FCT(2,:))'-num(1:28)']
ffT2 = toc;
ffT(n,1) = ffT2-ffT1;
    end
end

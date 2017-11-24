function rtTraCKerTrace
% minTraceLength = 2
% traceJmpForCombination = 1
cd('E:\MATLAB\TIRFcalibration\data\Ata01_5_125X100Y50x50_realtime')    
    %% files
    cd waSeq\tracker\
    cfg_ = load('..\..\cfgRT');
    cfg = cfg_.cfg;
    ndigit = cfg.ndigit; % # of digits for sequence number
    label = cfg.label; % # of digits for sequence number
    
    outDIR = 'rtData\';
    if exist(outDIR), rmdir(outDIR,'s'); end
    mkdir(outDIR)
    
    dn = cfg.sptReAppearTime;
    
    if dn>2, error('update the code (sptReAppearTime)'); end
    sptJmpForTracing = cfg.sptJmpForTracing;
    
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
    while (1)
        %% output filename
        digitTXT = eval(['sprintf(' digitFormat ',n)'] );
        traceDataFileNm = [outDIR 'traceData_' label '_' digitTXT '.mat'];
        
        %% load data
        while (1) % wait for update
            posFN = rdir(['posData-coeff*_' label '_' digitTXT '.mat']);
            if isempty(posFN)
                [fdbck] = funcFeedback(cfg.msgTXT,fdbck,fcall);
                if fdbck.inStop, break;  end % STOP
            else
                break; % continue
            end 
        end
        posfn = posFN.name;
        % scatter(XC,YC,'.')
        
        posfn_=load(posfn);
        XC = posfn_.X;
        YC = posfn_.Y;
        INTC = posfn_.INT;
        X = [X XC];
        Y = [Y YC];
        INT = [INT INTC];
        
        nSPOTs(n) = length(XC); % number of spots in that frame
        %ixSpt = ixSptFrm(n):ixSptFrm(n)+nSPOTs(n)-1; % indices of the spots in the current frame
        %frmNoSpot(ixSpt) = n;
        ixSptFrm(n+1) = ixSptFrm(n)+nSPOTs(n); % first spot in the frame
        if n == 1, n = 2; continue; end

        %% FIND OUT THE TRACES
        %%initialize
        % 1: 1st frame
        % 2: number of frames
        % 3: position in the trace array
        % 4-6: mean x, y , int
        % 7 : std deviation from the center

        for m = 1:nSPOTs(n) % number of spots
            Xc=XC(m);       % current
            Yc=YC(m);       % 
            n1 = n-dn;
            n2 = n-1;
            if n1<1, n1 = 1; end
            if n2<1, n2 = 1; end
            NS = fliplr(n1:n2); % search frames
            
            % check for addition to a prev trace
            isAdd2Trace = 0;
            for l=1:numel(NS) % earlier frames
                if isAdd2Trace, continue; end
                if isempty(trInf), continue; end
                ns = NS(l);
                last = trInf(:,1)+trInf(:,2)-1; % last frames of the traces
                ixf = find(last == ns); % traces in the current frame
                if isempty(ixf), continue; end
                tx = cellfun(@(v) v(end), TraceX(ixf));
                ty = cellfun(@(v) v(end), TraceY(ixf));
                dif=sqrt((Xc-tx).^2+(Yc-ty).^2);
                [~,ixmin] = min(dif);
                if dif(ixmin) < sptJmpForTracing % add to trace
                    tix = ixf(ixmin);
                    xm = []; ym = [];
                    if l > 1 % fill the jump
                        xm = (TraceX{tix}(end)+Xc)/2;
                        ym = (TraceY{tix}(end)+Yc)/2;
                    end
                    TraceX{tix} = [TraceX{tix} xm Xc];
                    TraceY{tix} = [TraceY{tix} ym Yc];
                    %TraceX{tix} = [TraceX{tix} nan(1,l-1) Xc];
                    %TraceY{tix} = [TraceY{tix} nan(1,l-1) Yc];
                    trInf(tix,2) = trInf(tix,2) + l;
                    isAdd2Trace = 1;
                end
            end

            % check for new trace
            isNewTrace = 0;
            ixSptTr = nan(1,dn+1);
            ixSptTr(end) = ixSptFrm(n)+m-1;
            if ~isAdd2Trace
                for l=1:numel(NS) % earlier frames
                    if isNewTrace, continue; end
                    ns = NS(l);
                    ixSptSearchFrm = ixSptFrm(ns):ixSptFrm(ns+1)-1;
                    dif=sqrt((Xc-X(ixSptSearchFrm)).^2+(Yc-Y(ixSptSearchFrm)).^2);
                    [~,ixmin] = min(dif);
                    if dif(ixmin) < sptJmpForTracing % add to trace
                        nb = n-ns; % backwards
                        ixSptTr(end-nb) = ixSptFrm(ns)+ixmin-1; % index of the spot in X array
                        isNewTrace = 1;
                    end
                end
                if isNewTrace
                    prevFrm = find(ixSptTr>0,1); % 
                    ixSptTr = ixSptTr(prevFrm:end);
                    xm = []; ym = [];
                    if sum(isnan(ixSptTr)) % fill the jump
                        xm = (X(ixSptTr(1))+X(ixSptTr(end)))/2;
                        ym = (Y(ixSptTr(1))+Y(ixSptTr(end)))/2;
                    end
                    TraceX{end+1} = [X(ixSptTr(1)) xm X(ixSptTr(end))];
                    TraceY{end+1} = [Y(ixSptTr(1)) ym Y(ixSptTr(end))];
                    %TraceX{end+1} = [X(ixSptTr(1)) nan(1,sum(isnan(ixSptTr))) X(ixSptTr(end))];
                    %TraceY{end+1} = [Y(ixSptTr(1)) nan(1,sum(isnan(ixSptTr))) Y(ixSptTr(end))];
                    trInf(end+1,1) = n-nb;
                    trInf(end,2) = 2;
                    X(ixSptTr(1)) = inf;
                    X(ixSptTr(end)) = inf;
                    Y(ixSptTr(1)) = inf;
                    Y(ixSptTr(end)) = inf;
                end
            end
        end
        cc= 3;
        n = n + 1;
        save(traceDataFileNm,'TraceX','TraceY','trInf')
        
        
        continue;

           

        
        %% PLOT1 : fill jumps
        if ~exist(traceJmplessDataFileNm,'file')

            TraceX2 = TraceX;
            TraceY2 = TraceY;        
            for tr = 1:size(trInf,1) % number of traces
                %% fill gaps (missing frames in jumpy traces)
                traceX = TraceX(trInf(tr,3):trInf(tr,3)+trInf(tr,2)-1);
                traceY = TraceY(trInf(tr,3):trInf(tr,3)+trInf(tr,2)-1);
                traceX(traceX<0) = 0.001; % no negative values
                traceY(traceY<0) = 0.001;

                ind = find(traceX~=0);
                frst = min(ind); last = max(ind);
                %iJump = find(trace(frst:last)==0)+frst-1;

                jmp = abs((traceX(frst:last)>0)-1);
                if sum(jmp > 0)
                    bnd = bwboundaries(jmp,'noholes');
                    for i = 1:numel(bnd) % for each jump
                        temp = bnd{i}+frst-1; % jump boundaries
                        jb = temp(:,2); clear temp;
                        mx= max(jb); mn=min(jb);
                        jL = mx-mn+2; % length

                        jSx = traceX(mx+1)-traceX(mn-1);% size
                        jSy = traceY(mx+1)-traceY(mn-1);% size
                        jsX = jSx/jL;% step
                        jsY = jSy/jL;% step
                        jVx = traceX(mn-1)+jsX*(1:jL-1);
                        jVy = traceY(mn-1)+jsY*(1:jL-1);
                        traceX(mn:mx)=jVx;
                        traceY(mn:mx)=jVy;
                        TraceX2(trInf(tr,3):trInf(tr,3)+trInf(tr,2)-1) = traceX;
                        TraceY2(trInf(tr,3):trInf(tr,3)+trInf(tr,2)-1) = traceY;
                    end
                end  
            end
            save(traceJmplessDataFileNm,'TraceX2','TraceY2')
        end
        
        
    end % while loop
    
    
end

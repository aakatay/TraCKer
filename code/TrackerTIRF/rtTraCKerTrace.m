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
            n1 = n-dn-1;
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
                    TraceX{tix} = [TraceX{tix} nan(1,l-1) Xc];
                    TraceY{tix} = [TraceY{tix} nan(1,l-1) Yc];
                    trInf(tix,2) = trInf(tix,2) + l;
                    isAdd2Trace = 1;
                end
            end

            % check for new trace
            isNewTrace = 0;
            ixSptTr = nan(1,dn+2);
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
                    TraceX{end+1} = [X(ixSptTr(1)) nan(1,sum(isnan(ixSptTr))) X(ixSptTr(end))];
                    TraceY{end+1} = [Y(ixSptTr(1)) nan(1,sum(isnan(ixSptTr))) Y(ixSptTr(end))];
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

           
        %% save results
        trInf = trInf(sum(trInf,2)>0,:);
%        weight = sum(TraceX>0,2);
%        nHist = hist(weight,max(weight)-min(weight));
%        histTraceLen  = [zeros(min(weight)-1) nHist];


        %%= 2 ==================================================================================================
        TraceX2 = TraceX;
        TraceY2 = TraceY;
        TraceINT2 = TraceINT;
        trInf2 = trInf;
        
        %% COMBINE TRACES
        Boy2 = size(trInf2,1);
        i=0;
        while i < size(trInf2,1)-1 % each trace
            i=i+1;
        %     TraceDif1=0;TraceDif2=0;TraceDif3=0;
            LastElement = trInf2(i,3)+trInf2(i,2)-1; % (position in the trace array)+(number of frames)
            LastBefore=LastElement-1;
            j=i;
            while j < size(trInf2,1)
                j=j+1;
                FirstElement=trInf2(j,3); % (position in the trace array)
                FirstAfter=FirstElement+1;
                diffTime = FirstElement-LastElement; % [frames]
                if (0 <= diffTime)  && (diffTime <= 2) 
                    TraceDif1=sqrt([TraceX2(LastElement)-TraceX2(FirstElement)]^2+[TraceY2(LastElement)-TraceY2(FirstElement)]^2);
                    TraceDif2=sqrt([2*TraceX2(LastElement)-TraceX2(LastBefore)-TraceX2(FirstElement)]^2+[2*TraceY2(LastElement)-TraceY2(LastBefore)-TraceY2(FirstElement)]^2);
                    TraceDif3=sqrt([TraceX2(LastElement)-2*TraceX2(FirstElement)+TraceX2(FirstAfter)]^2+[TraceY2(LastElement)-2*TraceY2(FirstElement)+TraceY2(FirstAfter)]^2);
                    TraceDif=[TraceDif1,TraceDif2,TraceDif3];
                    %TraceDif=sqrt((TraceX2(i,LastElement)-TraceX2(j,FirstElement))^2+(TraceY2(i,LastElement)-TraceY2(j,FirstElement))^2);
                    if min(TraceDif)<traceJmpForCombination % combine traces
                        isOverlapping = 0;
                        isSpace = 0;
                        if diffTime == 2 % overlapping
                            isSpace = 1; % one frame is empty 
                        elseif diffTime == 0 % overlapping
                            isOverlapping = 1;
                            k1 = LastElement;
                            k2 = FirstElement;
                            TraceX2(k1)= [TraceX2(k1)*TraceINT2(k1)+TraceX2(k2)*TraceINT2(k2)]/[TraceINT2(k1)+TraceINT2(k2)];
                            TraceY2(k1)= [TraceY2(k1)*TraceINT2(k1)+TraceY2(k2)*TraceINT2(k2)]/[TraceINT2(k1)+TraceINT2(k2)];
                            TraceINT2(k1)= [TraceINT2(k1)*TraceINT2(k1)+TraceINT2(k2)*TraceINT2(k2)]/[TraceINT2(k1)+TraceINT2(k2)];
                        end
                        trInf2(i,2) = trInf2(i,2)+trInf2(j,2)-isOverlapping+isSpace; % number of frames
                        % join the traces
                        fr1   = LastElement+1+isSpace; % extension of the first trace
                        fr2   = LastElement+isSpace+trInf2(i,2);
                        fr1_2 = FirstElement+isOverlapping; % coord. of the second trace
                        fr2_2 = FirstElement+trInf2(j,2)-1;
                        TraceX2(fr1:end) = [TraceX2(fr1_2:fr2_2) TraceX2(fr1:fr1_2-1) TraceX2(fr2_2+1:end) ];

                        LastElement=trInf2(i,1)+trInf2(i,2)-1;
                        LastBefore=LastElement-1;
                        trInf2(j:end-1,:) = trInf2(j+1:end,:);
                        trInf2 = trInf2(1:end-1,:);

                        % update trace info of the following traces
                        trInf2(j:end,3)=trInf2(j:end,3)-isOverlapping+isSpace;
                        trInf2(i+1:j-1,3)=trInf2(i+1:j-1,3) + fr2_2-fr1_2+isSpace;
                        % update trace values
                        com=com+1;
                    end
                end
            end % i <= Boy2-1
        end
        save(traceDataFileNm,'TraceX2','TraceY2','TraceINT2','trInf2','frmNoTrace')
        
        n = n + 1;
    end
end

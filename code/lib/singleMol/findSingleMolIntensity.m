function [intSingleMolperFrame, intSingleMolperTrace] = findSingleMolIntensity(traceData0FN,acqTime)

        
    maxTrLen=10;
    frBsz = 2000; % frame block size
    bw=2000; % intensity bin width
    minTimeReq = 360; % [sec] (uses late frames
    minTimeReq = 120; % [sec] (uses late frames)
    isSingleFrame = 1; % ow. displays trace total intensity
    isDisp = 1;
    
    load(traceData0FN,'TraceINT','frmNoTrace','trInf');
    %remove long traces
    nTr0 = size(trInf,1);
    trInf = trInf(trInf(:,2)<=maxTrLen,:);
    nTr1 = size(trInf,1);
    disp(sprintf('%.02f%% of the traces are long and excluded in histogram analysis',(nTr0-nTr1)/nTr0*100))

    % k=1: % int per frame
    % k=2: % total trace intensity
    for k = 1:2 
    %for k = 1
        if k==1
            INT = TraceINT; % int per frame
            frm = frmNoTrace;
        else % trace
            frmLen = trInf(:,2); % trace length [frames]
            INT = trInf(:,6).*frmLen; % total trace intensity
            frm = trInf(:,1); % 1st frame trace starts
        end
        
        Nfr = max(frm);
        nB = floor(Nfr/frBsz);
        Nfr2 = nB*frBsz;
        if acqTime*Nfr2<minTimeReq, warning('short acquisition duration to detect single molecule intensity (min:6min)'); end;
        [N,xed] = histcounts(INT,'Binwidth',bw);
        ixmxALL = floor(max(xed)/bw);
        ixmnALL = floor(min(xed)/bw)+1;
        xc = xed(1:end-1)+bw/2; % centers
        xedALL = xed;
        TraceIntHist = zeros(nB,ixmxALL);
        j=1;
        for i = 0:frBsz:Nfr
            i1 = i;
            i2 = i+frBsz;
            ix = find(((i1<frm).*(frm<i2)));
            [N,xed] = histcounts(INT(ix),'Binwidth',bw);
            %plot(N)
            ixmx = floor(max(xed)/bw);
            ixmn = floor(min(xed)/bw)+1;
            TraceIntHist(j,ixmn:ixmx)=N;
            j=j+1;
        end

        % display
        % 3D histogram image 
        %TraceIntHist(:,1)=0;
        TraceIntHistNorm=repmat(1./sum(TraceIntHist,2),1,size(TraceIntHist,2)).*TraceIntHist;
        if isDisp
            figure(100);
            imagesc(TraceIntHistNorm);
            xedDisp = 0:bw:max(xedALL);
            set(gca,'XTick',0.5:numel(xedDisp)+0.5,'XTickLabel',xedDisp,'XTickLabelRotation',90)
            set(gca,'YTick',0.5:Nfr+0.5,'YTickLabel',0:frBsz:Nfr+frBsz)
            xlabel('intensity'); ylabel('frames')
            title('single molecule intensity histogram by time')
        end

        % max bin intensity plot
        TraceIntHistNormSmth = conv2(TraceIntHistNorm,[1 1 1]/3,'same');
        [~,intvsframe] = max(TraceIntHistNormSmth,[],2);
        %intvsframe= intvsframe+ixmnALL;
        y = intvsframe*bw-bw/2;
        x = (1:nB+1)*frBsz;

        % determine single mol intensity
        ix = find(x*acqTime>minTimeReq-60);
        TraceIntHistLast = sum(TraceIntHistNormSmth(ix(1):end,:),1);
        [~,ixmax] = max(TraceIntHistLast,[],2);
        intSingleMol(k) = ixmax*bw-bw/2;

        if isDisp
            figure(101);
            f=fit(x',y,'poly2');
            plot(f,x',y,'o')
            xlabel('frames')
            ylabel('single mol. intensity')
            grid minor
            line([x(ix(1)) x(ix(1))],[min(y) max(y)])
            imline(gca,[0 max(x)],[intSingleMol(k) intSingleMol(k) ])

            % 2D histogram
            %figure(102)
            %bar(sum(TraceIntHist,1))
            %grid minor
        end
    end
    intSingleMolperFrame = intSingleMol(1);
    intSingleMolperTrace = intSingleMol(2);

end
clear all;
close all;

sptReAppearTime = 2;
sptJmpForTracing = 1;
minTraceLength = 2;

        % input file
        xyzDataGausFileNm = 'xyzDataGaus-coeff4070_002.mat';
        load(xyzDataGausFileNm); 
        % output files
        traceDataFileNm0 = 'traceData0-coeff4070_002.mat';
        traceJmplessDataFileNm = 'traceJmplessData-coeff4070_002.mat';
        
        
        %% fill zero values
        [iy, ix] = find(X>0); %debugSpt = [X(X>0) Y(Y>0) ix];
        while ~isempty(find(ixSptFrm==0))
            ixSptFrm(find(ixSptFrm==0))=ixSptFrm(find(ixSptFrm==0)+1)
        end
                                
        
        %% initialize
        isCropFrames = 0;
        if isCropFrames
            f1=1;f2=300;
            X = X(:,f1:f2);
            Y = Y(:,f1:f2);
        end
        Xilk=X;
        Yilk=Y;
        nSpots = numel(X); % number of spots
        Frames = numel(ixSptFrm)-1; % number of frames
        nTRg = round(nSpots/3); % initial guess for number of traces
        trInf = zeros(nTRg,3); % trace info, memory allocation
        % 1: 1st frame
        % 2: number of frames
        % 3: position in the trace array
        % 4-6: mean x, y , int
        % 7 : std deviation from the center

        %% LOOP
        p=1;
        f=0;
        tic
        isDebug = 0;
        ixTr = 1;
        h = waitbar(0,'Finding the traces...');
        for k=1:Frames-1 % number of frames
            if isDebug, disp(sprintf('=frame#:%i ',k)); end
            Boy = ixSptFrm(k+1) - ixSptFrm(k); % num spots in the frame
            for m_=1:Boy % number of spots
                clear ixSptTr;
%                dif=Inf(Boy,Frames-k+1);
%                difbin=zeros(Boy,Frames-k+1);
                m = m_ + ixSptFrm(k) -1;
                ixSptTr = m;    % first element of the trace
                AslX=X(m);       % last X value in the tracking
                AslY=Y(m);       % last Y value in the tracking
                if AslX == inf
                    continue;
                end
                if isDebug, disp(sprintf('==spot#:%i ',m)); end
                quit = 0;
                ll=k;
                isNewTrace = 0;
                for l=k+1:Frames % later frames
                    if isDebug, disp(sprintf('===check frame#:%i',l)); end
                    if quit
                        break; 
                    end
                    ixSptSearchFrm = ixSptFrm(l):ixSptFrm(l+1)-1;
                    dif=sqrt((AslX-X(ixSptSearchFrm)).^2+(AslY-Y(ixSptSearchFrm)).^2);
                    [v,n] = min(dif);
                    BOY = Boy;
                    BOY = 1;
                    %for n=1:BOY % all spots
                    if l-ll>sptReAppearTime % if the next frame where a spots re-appears in sptJmp distance is 4 frames apart ignores it.
                        %if sum(sum(difbin(:,l-2:l-1))) == 0, 
                            quit =1;
                            break, 
                        %end
                    end

                    if isDebug, disp(sprintf('===check spot#:%i ',n)); end
                    if dif(n) < sptJmpForTracing; % add to trace
                        %if dif(n,l)==min(dif(:,l))     
                            ll=l;
                            ixSptTr(l-k+1) = ixSptFrm(l)+n-1; % index of the spot in X array
                            AslX=Xilk(ixSptTr(l-k+1));
                            AslY=Yilk(ixSptTr(l-k+1));
                            isNewTrace = 1;
                            if isDebug, disp(sprintf('===add__ spot#:%i, #%i ',n,l-k+1)); end
                            %break
                        %end
                    end
                    %end % spots
                end % frames
                
                if isNewTrace
                    clear tracex tracey traceint;
                    tracex(ixSptTr>0) = Xilk(ixSptTr(ixSptTr>0));
                    tracey(ixSptTr>0) = Yilk(ixSptTr(ixSptTr>0));
                    traceint(ixSptTr>0) = INT(ixSptTr(ixSptTr>0));
                    X(ixSptTr(ixSptTr>0))=Inf;
                    Y(ixSptTr(ixSptTr>0))=Inf;
                    num=numel(find(tracex>0)); % trace length
                    if num>=minTraceLength % if # of data points larger than minTraceLength, than saves as a trace
                        pos=find(tracex>0);
                        ilk=zeros(1,num+1);
                        son=zeros(1,num+1);
                        ilk(1:1:num)=pos(1:1:num);
                        son(2:1:num+1)=pos(1:1:num);
                        fark=ilk-son;
                         %if numel(find(fark==1))>2 % # of consecutive data points
                        p2 = p+numel(tracex)-1; % position of the spot in trace array
                        TraceX(p:p2)= tracex;
                        TraceY(p:p2)= tracey;
                        TraceINT(p:p2)= traceint;   
                        frmNoTrace(p:p2) = k:k+numel(tracex)-1; %% frame numbers of each spot in the trace
                        trInf(ixTr,1) = k; % first frame of the trace
                        trInf(ixTr,2) = numel(tracex); % number of frames (length)
                        trInf(ixTr,3) = p; % position in the trace arrays
                        trInf(ixTr,4) = mean(tracex(tracex>0));
                        trInf(ixTr,5) = mean(tracey(tracex>0));
                        trInf(ixTr,6) = mean(traceint(tracex>0));
                        trInf(ixTr,7) = 0; % std(sqrt((tracex-trInf(ixTr,4)).^2+(tracey-trInf(ixTr,5)).^2)); % RMS deviation
                        trInf(ixTr,8) = findMaxDist(tracex,tracey); % MAX distance
                        %trInf(ixTr,8) = gQF; % gaussian quality factor

                        p=p2+1;
                        ixTr = ixTr + 1; % trace index
                         %end
                    end
                end


            end
           waitbar(k / Frames)
           if exist('_stopRunning-ON','file')
               break
           end
        end
        trInf = trInf(sum(trInf,2)>0,:);
        TrackTime = toc;
%        weight = sum(TraceX>0,2);
%        nHist = hist(weight,max(weight)-min(weight));
%        histTraceLen  = [zeros(min(weight)-1) nHist];
        
        save('TrackTime','TrackTime');
        close(h)
        %return; 
        save(traceDataFileNm0,'TraceX','TraceY','TraceINT','trInf','frmNoTrace'); 

        
        %% combine
        
    if 0 && ~exist(traceDataFileNm) && isCombTraces % speed and trace combination
        %% COMBINE TRACES
        load(traceDataFileNm0)
        [m n]=size(TraceX);
        h = waitbar(0,'Combining the traces...');
        disp('Combining the traces...');
        load('xyzDataGaus-coeff941_000.mat', 'ixSptFrm')
        Boy2 = size(trInf,1);
        Frames = numel(ixSptFrm)-1;
        tic;
        com=0;
        i=0;
        while i < size(trInf,1)-1 % each trace
            i=i+1;
        %     TraceDif1=0;TraceDif2=0;TraceDif3=0;
            LastElement = trInf(i,3)+trInf(i,2)-1; % (position in the trace array)+(number of frames)
            LastBefore=LastElement-1;
            j=i;
            while j < size(trInf,1)
                j=j+1;
                FirstElement=trInf(j,3); % (position in the trace array)
                FirstAfter=FirstElement+1;
                diffTime = FirstElement-LastElement; % [frames]
                if (0 <= diffTime)  && (diffTime <= 2) 
                    TraceDif1=sqrt([TraceX(LastElement)-TraceX(FirstElement)]^2+[TraceY(LastElement)-TraceY(FirstElement)]^2);
                    TraceDif2=sqrt([2*TraceX(LastElement)-TraceX(LastBefore)-TraceX(FirstElement)]^2+[2*TraceY(LastElement)-TraceY(LastBefore)-TraceY(FirstElement)]^2);
                    TraceDif3=sqrt([TraceX(LastElement)-2*TraceX(FirstElement)+TraceX(FirstAfter)]^2+[TraceY(LastElement)-2*TraceY(FirstElement)+TraceY(FirstAfter)]^2);
                    TraceDif=[TraceDif1,TraceDif2,TraceDif3];
                    %TraceDif=sqrt((TraceX(i,LastElement)-TraceX(j,FirstElement))^2+(TraceY(i,LastElement)-TraceY(j,FirstElement))^2);
                    if min(TraceDif)<traceJmpForCombination % combine traces
                        isOverlapping = 0;
                        isSpace = 0;
                        if diffTime == 2 % overlapping
                            isSpace = 1; % one frame is empty 
                        elseif diffTime == 0 % overlapping
                            isOverlapping = 1;
                            k1 = LastElement;
                            k2 = FirstElement;
                            TraceX(k1)= [TraceX(k1)*TraceINT(k1)+TraceX(k2)*TraceINT(k2)]/[TraceINT(k1)+TraceINT(k2)];
                            TraceY(k1)= [TraceY(k1)*TraceINT(k1)+TraceY(k2)*TraceINT(k2)]/[TraceINT(k1)+TraceINT(k2)];
                            TraceZ(k1)= [TraceZ(k1)*TraceINT(k1)+TraceZ(k2)*TraceINT(k2)]/[TraceINT(k1)+TraceINT(k2)];
                            TraceINT(k1)= [TraceINT(k1)*TraceINT(k1)+TraceINT(k2)*TraceINT(k2)]/[TraceINT(k1)+TraceINT(k2)];
                        end
                        trInf(i,2) = trInf(i,2)+trInf(j,2)-isOverlapping+isSpace; % number of frames
                        % join the traces
                        fr1   = LastElement+1+isSpace; % extension of the first trace
                        fr2   = LastElement+isSpace+trInf(i,2);
                        fr1_2 = FirstElement+isOverlapping; % coord. of the second trace
                        fr2_2 = FirstElement+trInf(j,2)-1;
                        TraceX(fr1:end) = [TraceX(fr1_2:fr2_2) TraceX(fr1:fr1_2-1) TraceX(fr2_2+1:end) ];

                        %frmNoTrace(fr1:end) = 
                        LastElement=trInf(i,1)+trInf(i,2)-1;
                        LastBefore=LastElement-1;
                        trInf(j:end-1,:) = trInf(j+1:end,:);
                        trInf = trInf(1:end-1,:);

                        % update trace info of the following traces
                        trInf(j:end,3)=trInf(j:end,3)-isOverlapping+isSpace;
                        trInf(i+1:j-1,3)=trInf(i+1:j-1,3) + fr2_2-fr1_2+isSpace;
                        % update trace values

                        %trInf(i,4) = mean(tracex(tracex>0));
                        %trInf(i,5) = mean(tracey(tracex>0));
                        %trInf(i,6) = mean(traceint(tracex>0));
                        %trInf(i,7) = std(sqrt((tracex-trInf(ixTr,4)).^2+(tracey-trInf(ixTr,5)).^2)); % RMS deviation                            

                        % update 

                        com=com+1;
                    end
                end
            end % i <= Boy2-1
        end

        save(traceDataFileNm,'TraceX','TraceY','TraceZ','TraceINT','TraceSpeed','trInf','frmNoTrace')
        delete(traceJmplessDataFileNm)
    end
    

    %return
    %% PLOT1 : fill jumps
    if ~exist(traceJmplessDataFileNm,'file')
        disp('preparing plot data...')
        
        % trInf :
        % 1: first frame of the trace
        % 2: number of frames (length)
        % 3: position in the trace arrays
        
        
        tic
        hWBfillgaps = waitbar(0,'filling gaps');
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
            waitbar(tr/size(TraceX,1))
        end
        close(hWBfillgaps)
        toc
        save(traceJmplessDataFileNm,'TraceX2','TraceY2')
    end

%return % un comment to stop


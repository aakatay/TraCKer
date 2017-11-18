function [ traceVol ] = dispTraceVol( traceDataPath(i),rotAngle(i),tiltAngle(:,i) )
%DISPTRACEVOL displays the 3D time lag evolution of the spots
%   imgStcks : 4D data of 3D time lag data
    %% read trace coord.
    load 
    traceCoor(1) = TraceX;
    traceCoor(2) = TraceY;
    traceCoor(3) = TraceZ;
    traceCoorDS = mapTraceCoor(traceCoor,);
    [uv sv] = memory;
    
    memAv = sv.PhysicalMemory.Available; % available memory in bytes
    memAv/1024/1024/1024 % available memory in Gbytes
    
   % determine the max.number of frames that canbe processed without paging
   % blocks to be processed separately
   fileSize = imageInfo(1).Width*imageInfo(1).Height*2; % bytes
   maxNfrm = floor(memAv/fileSize);
   nBlcks = ceil(frmN/maxNfrm)
   blckNfrm = ceil(nFrm/nBlcks); % # of 
   = 
   frm
    if exist('stack.mat')
        load stack;
    else
        for k=1:Frames*ratFrm
            for l=1:StackNum

                DigitDiff=floor(log10(Frames))-floor(log10(k));   
                %fprintf('Frame:%d Stack:%d \n',k,l)

                if k == 1
                 DigitDiff=floor(log10(Frames));
                end

                if DigitDiff == 0
                Stack(:,:,l,k)=uint16(imread(['stack_' int2str(k) '.tif'],l));
                end

                if DigitDiff == 1
                Stack(:,:,l,k)=uint16(imread(['stack_0' int2str(k) '.tif'],l));
                end

                if DigitDiff == 2
                Stack(:,:,l,k)=uint16(imread(['stack_00' int2str(k) '.tif'],l));
                end

                if DigitDiff == 3
                Stack(:,:,l,k)=uint16(imread(['stack_000' int2str(k) '.tif'],l));
                end
            end
            h = waitbar(k/Frames/ratFrm);
        end
    end    

end


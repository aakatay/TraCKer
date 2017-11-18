% converts the linear indexed data format to array indexed format
function [TraceX TraceY TraceInt trInf frmNoTrace] = convertArray2Linear(arrayData)

nTrace = size(arrayData,1);  % number of traces
p=1;
for i = 1:nTrace
    frm = arrayData(i,arrayData(i,:,1)>0,1);
    % correct if frame numbers are increasing
    [~, ix]=max(frm);frm = frm(1:ix);   
    
    nfrm = frm(end)-frm(1)+1;
    nElem = numel(frm); % number of elements
    nzElem = 1:nElem;
    traceX = zeros(nfrm,1);
    traceY = traceX;
    traceInt = traceY;
    
    arrayData(i,2);
    traceX(frm-frm(1)+1) = arrayData(i,nzElem,2);
    traceY(frm-frm(1)+1) = arrayData(i,nzElem,3);
    traceInt(frm-frm(1)+1) = arrayData(i,nzElem, 4);
    
    ind = find(traceX~=0);
    frst = min(ind); last = max(ind);
    
    jmp = abs((traceX'>0)-1);
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
        end
    end  
    
    
    trInf(i,1) = frm(1);
    trInf(i,2) = nfrm;
    trInf(i,3) = p;
    trInf(i,4) = mean(traceX);
    trInf(i,5) = mean(traceY);
    trInf(i,6) = mean(traceInt);
    trInf(i,7) = std(sqrt((traceX-trInf(i,4)).^2+(traceY-trInf(i,5)).^2));
    trInf(i,8) = findMaxDist(traceX',traceY'); 
    
    frmNoTrace(p:p+nfrm-1) = frm(1):frm(1)+nfrm-1;
    TraceX(p:p+nfrm-1)= traceX;
    TraceY(p:p+nfrm-1)= traceY;
    TraceInt(p:p+nfrm-1)= traceInt;
    p = p + nfrm;
        % 1: 1st frame
        % 2: number of frames
        % 3: position in the trace array
        % 4-6: mean x, y , int
        % 7 : std deviation from the center
end
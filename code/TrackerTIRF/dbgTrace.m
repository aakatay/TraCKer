% display traces with TrackerTIRF
clear all;
close all;




tx = 1; % trace indices
N = 100;
for i = 1 : N 
    tx = i;
    % (2) TrackerTIRF.m
    load traceDbg0; % trInf [traceData0*.mat]
    trInf2 = trInf; clear trInf;
    load traceDbg; % TraceX2 TraceY2 (traceJmplessData*.mat)

    ix1 = trInf2(tx,3);
    ix2 = trInf2(tx,2)+ix1-1;

    x2 = TraceX2(ix1:ix2)';
    y2 = TraceY2(ix1:ix2)';

    % (1) rtTraCKerTrace.m
    load('traceData_shft3.5_0017.mat') 
    x = TraceX{tx}';
    y = TraceY{tx}';

    n2 = numel(x2);
    n1 = numel(x);
    if n1 ~= n2
        nmx = max(n1,n2);
        a = nan(nmx,4);
        a(1:n2,1)=x2;
        a(1:n1,2)=x;
        a(1:n2,3)=y2;
        a(1:n1,4)=y
        cccc= 3;
    else
        [x2 x y2 y]
    end
        
    
    
end
% display traces with TrackerTIRF
clear all;
close all;
isConverted =1;
    % (1) rtTraCKerTrace.m
    if isConverted % use converted data
        load('traceDataRT.mat')
    else
        load('traceData_shft3.5_0017.mat') 
    end
    TraceXrt = TraceX; TraceYrt = TraceY; trInfRT = trInf;
    
    % (2) TrackerTIRF.m (from C:\MATLAB\TraCKer\code\TrackerTIRF)
    load traceDbg0; % trInf [traceData0*.mat]
    load traceDbg; % TraceX2 TraceY2 (traceJmplessData*.mat)



tx = 1; % trace indices
N = 100;
for i = 1 : N 
    tx = i; % trace index

    % (1) rtTraCKerTrace.m
    if isConverted
        ix1 = trInfRT(tx,3);
        ix2 = trInfRT(tx,2)+ix1-1;

        x = TraceXrt(ix1:ix2)';
        y = TraceYrt(ix1:ix2)';
    else
        x = TraceXrt{tx}';
        y = TraceYrt{tx}';
    end

    % (2) TrackerTIRF.m
    ix1 = trInf(tx,3);
    ix2 = trInf(tx,2)+ix1-1;

    x2 = TraceX2(ix1:ix2)';
    y2 = TraceY2(ix1:ix2)';
        

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
        pause
    else
        [x2 x y2 y]
    end
        
    
    
end
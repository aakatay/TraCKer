% converts realtime data to bulk process data 
% + for TraCKerTIRFdisp.m
clear all;
close all;
fnTraceDataOUT = '..\..\..\traceDataRT';
%% load traces
    FN = dir('traceData_*.mat');
    n=  numel(FN);
    IX = [];

    for i = 1:n
        fn = FN(i).name;
        IX  = [IX; str2num(fn(end-7:end-4))];
    end
    [~,IXs] =sort(IX);
    FN = FN(IXs);

    for i = 1:n % check 
        fn = FN(i).name;
        IX  = [IX; str2num(fn(end-7:end-4))];
    end
    
    load(FN(end).name)
    tx = TraceX; ty = TraceY; tr = trInf;
    clear TraceX TraceY;
    nt = size(trInf,1);
    ix = 1;
    for i = 1:nt % each trace
        x = tx{i};
        y = ty{i};
        L = numel(x); % number of frames
        TraceX(ix:ix+L-1) = x;
        TraceY(ix:ix+L-1) = y;
        trInf(i,3) = ix;
        ix = ix + L;
    end
    
    save(fnTraceDataOUT,'TraceX','TraceY','trInf');
    cd('..\..\..\')
if exist('fileTraceData')
    save(fileTraceData,'TraceX','TraceY','TraceZ','TraceINT');
else
    save('traceData.mat','TraceX','TraceY','TraceZ','TraceINT');
end
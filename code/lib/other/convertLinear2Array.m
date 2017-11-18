% converts the linear indexed data format to array indexed format
function [TraceXarray TraceYarray] = convertLinear2Array(TraceX,TraceY,trInf)

nTrace = size(trInf,1);  % number of traces
for i = 1:nTrace
    TraceXarray(i,trInf(i,1):trInf(i,1)+trInf(i,2)-1) = TraceX((trInf(i,3):trInf(i,3)+trInf(i,2)-1));
    TraceYarray(i,trInf(i,1):trInf(i,1)+trInf(i,2)-1) = TraceY((trInf(i,3):trInf(i,3)+trInf(i,2)-1));
end
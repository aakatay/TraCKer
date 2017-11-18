% converts the linear indexed data format to array indexed format
function [TraceXarray TraceYarray] = convertLinear2Array(X,Y,INT,ixSptFrm,frmNoSpot)

frames = size(ixSptFrm,1)-1;  % number of traces
for i = 1:frames
    Xarray(i,ixSptFrm(i,1):ixSptFrm(i+1,1)-1) = X((ixSptFrm(i,3):ixSptFrm(i+1,3)-1));
    Yarray(i,ixSptFrm(i,1):ixSptFrm(i+1,1)-1) = Y((ixSptFrm(i,3):ixSptFrm(i+1,3)-1));
    INTarray(i,ixSptFrm(i,1):ixSptFrm(i+1,1)-1) = INT((ixSptFrm(i,3):ixSptFrm(i+1,3)-1));
end
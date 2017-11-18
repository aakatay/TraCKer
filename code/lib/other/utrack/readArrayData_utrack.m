load('ProcessedTracks');


nTrace = size(tracks,2);
for i = 1:nTrace
    nEl(i) = numel(tracks(i).f);
end
mxEl = max(nEl);
Threshfxyc = zeros(nTrace,mxEl,4);
for i = 1:nTrace
    trLen = numel(tracks(i).f);
    Threshfxyc(i,1:trLen,1) = tracks(i).f;
    Threshfxyc(i,1:trLen,2) = tracks(i).x;
    Threshfxyc(i,1:trLen,3) = tracks(i).y;
    Threshfxyc(i,1:trLen,4) = tracks(i).A;
end
%fnameData = dir('acq*.mat');
%load(fnameData(1).name); % Threshfxyc, fxyc_struct

%Threshfxyc2 = permute(Threshfxyc,[3 1 2]);
%Threshfxyc2 = Threshfxyc2(:,:,[1 2 3 5] );
Threshfxyc2 = Threshfxyc;

traceLen = sum(im2bw(Threshfxyc2(:,:,1),0),2);
Threshfxyc3 = Threshfxyc2(traceLen>1,:,:); % discard 1 element traces

arrayData=Threshfxyc3;

[TraceX,TraceY,TraceInt,trInf,frmNoTrace] = convertArray2Linear(arrayData);

save('linearData','TraceX','TraceY','TraceInt','trInf','frmNoTrace');
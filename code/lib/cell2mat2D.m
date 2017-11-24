load('traceData_shft3.5_0017.mat')
nf = 17; % frames
nt = numel(TraceX); % traces

tracex = zeros(nt,nf);

for i = 1:nt
    tx = TraceX{i};
    tracex(i,1:numel(tx)) = tx;
    
end

tracex

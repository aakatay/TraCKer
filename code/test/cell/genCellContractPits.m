clear TraceX TraceY TraceInt
T = 20; % # frames
N = 128; % # px
TraceY = round(rand(N,1)).*round(rand(N,1)*N);
TraceY = repmat(TraceY,1,T);
TraceX(1,:) = 1:20;
TraceX = repmat(TraceX,N,1);

zeroIx = find(TraceY == 0);
TraceX(zeroIx) = 0;
TraceINT = ones(size(TraceX));
TraceINT(zeroIx) = 0;

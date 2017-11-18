% find the the amount of approach

%% find traces in prox
%function TrProxIx = findTracesInProx(TraceX,TraceY,T,x,y,BinSz)
% find traces positions closest to 
N = 10; % # of traces to calculate distance
Nsel = 5; % # of selected traces to evaluate in contraction 
Xc = (i-0.5)*BinSz;
Yc = (j-0.5)*BinSz;

% find continous traces
nonzeroT1 = TraceX(:,t)~=0;
nonzeroT2 = TraceX(:,t+1)~=0;
Tr = find((nonzeroT1+nonzeroT2) == 2); % both nonzero

% find neighbouring traces
dX = abs(TraceX(Tr,t)-Xc);
dY = abs(TraceY(Tr,t)-Yc);
[dP, trIx] = sort(dX+dY); % find distance parameter
p = 1:N; % trace index vector
dR1 = sqrt(   (TraceX(trIx(p)) - Xc).^2  +  (TraceY(trIx(p)) - Yc).^2   );
[dRs dRix] = sort(dR1);
TrProxIx = Tr(trIx(dRix(1:Nsel))); % indexes of the selected traces

% find approach param.
dXs = TraceX(TrProxIx,t)-Xc;
dYs = TraceY(TrProxIx,t)-Yc;
dRs = sqrt(dXs.^2 + dYs.^2);
dispVec = [ TraceX(TrProxIx,t+1) - TraceX(TrProxIx,t)   TraceY(TrProxIx,t+1) - TraceY(TrProxIx,t)    ]; % displacement vector
approachParam = sum(dispVec .* [dXs dYs],2) ./ sqrt(dXs.^2+dYs.^2); % dot product


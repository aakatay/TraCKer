% find the the amount of approach
%IN: TraceX, TraceY, TraceZ, Xc, Yc, Zc, pxszZvsXY, evalWinR, df

%% find traces in prox
%function TrProxIx = findTracesInProx(TraceX,TraceY,Xc,Yc)
% find traces positions closest to 
Tr = find(TraceX(:,f)>0);
% find neighbouring traces
dX = abs(TraceX(Tr,f)-Xc);
dY = abs(TraceY(Tr,f)-Yc);
dZ = abs(TraceZ(Tr,f)-Zc);
dR2 = sqrt(dX.^2+dY.^2+pxszZvsXY.^2*dZ.^2);
trIx = find(dR2<evalWinR); 
% indexes of the selected traces
ixTemp = trIx(find(dR2(trIx) ~= 0)); % remove self
TrProxIx = Tr(ixTemp);
dR = dR2(ixTemp);

% remove traces not continous in the following df frames
ixTemp = find(TraceX(TrProxIx,f+df)~=0);
TrProxIx2 = TrProxIx(ixTemp);
dR = dR(ixTemp);


    

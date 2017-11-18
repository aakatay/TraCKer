% 

clear

load('traceData.mat') % 'TraceX','TraceY','TraceZ','TraceINT','TraceT0'
% remove nan data
[Frm Trc]=find(~isnan(TraceX'));
yValues = TraceX(sub2ind(size(TraceX),Trc,Frm));
xValues = TraceY(sub2ind(size(TraceX),Trc,Frm));
zValues = TraceZ(sub2ind(size(TraceX),Trc,Frm));

scatter(xValues,yValues,zValues);
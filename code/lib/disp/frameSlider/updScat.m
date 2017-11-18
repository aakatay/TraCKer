if q, return; end;
    if isTrace == 1 % spots
        frmNo = frmNoSpot;
    elseif isTrace > 1 % traces
        frmNo = frmNoTrace;
    end
    
    CM = genColorMap('jet',numFrames);    
    ixSel = find((frmNo>=fr1) .* (frmNo<=fr2));     % events 
    ixSelTr = find((trInf2(:,1)>=fr1) .* ( (trInf2(:,1)+trInf2(:,2)) <=fr2));   % traces
    switch isTrace
        case 1 % spots
            xx = X(ixSel); yy = Y(ixSel);
            fr = frmNo(ixSel);
        case 2 % traces
            xx = TraceX(ixSel); yy = TraceY(ixSel);
            fr = frmNo(ixSel);
        case 3 % trace average
            xx = trInf2(ixSelTr,4); yy = trInf2(ixSelTr,5);
            fr = trInf2(ixSelTr,1);
        case 4 % trace recruit average (localized traces)
            ixSelTr2 = ixSelTr(trInf2(ixSelTr,6)>0.2);
            xx = trInf2(ixSelTr2,4); yy = trInf2(ixSelTr2,5);
            fr = trInf2(ixSelTr2,1);
    end
    fr = fr - min(fr) +1;
    %xx_=xx_(xx_>0); yy_=yy_(yy_>0);
     set(hScat,'CData',CM(fr,:),'XData',xx,'YData',yy);
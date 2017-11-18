function [funText] = genFitFunc(nspots)
% nspots : number of neighbouring spots
% ie if 8 spots :  1: bckgrnd 2: center spot 3-8: neigh. spots = > nspots=6

    nSpFitNonRegMax = 6; % 6 spots (neighbouring spots)
    nSpFitNonRegMax = nspots; % 6 spots (neighbouring spots)
    txFitSpNonRegFunc = 'funSpNonReg = @(c,x) c(1)'; % fit to non registered spots
    st = 0;
    for ixSpFit = 1:nSpFitNonRegMax
        txFitSpNonRegFunc = [txFitSpNonRegFunc sprintf('+c(%i)*exp(-intP*((x(:,1)-c(%i))/c(%i)/sqrt(2)).^2-intP*((x(:,2)-c(%i))/c(%i)/sqrt(2)).^2)',2+st,3+st,4+st,5+st,6+st)];
        st = st + 5;
    end
    funText = txFitSpNonRegFunc;
    return
end
    
        %% e.g.
% funSpNonReg = @(c,x) c(1)+ 26
% c(2)*exp(-intP*((x(:,1)-c(3))/c(4)/sqrt(2)).^2-intP*((x(:,2)-c(5))/c(6)/sqrt(2)).^2)+ 85
% c(7)*exp(-intP*((x(:,1)-c(8))/c(9)/sqrt(2)).^2-intP*((x(:,2)-c(10))/c(11)/sqrt(2)).^2)+ 87
% c(12)*exp(-intP*((x(:,1)-c(13))/c(14)/sqrt(2)).^2-intP*((x(:,2)-c(15))/c(16)/sqrt(2)).^2)+ 90
% c(17)*exp(-intP*((x(:,1)-c(18))/c(19)/sqrt(2)).^2-intP*((x(:,2)-c(20))/c(21)/sqrt(2)).^2)+ 90
% c(22)*exp(-intP*((x(:,1)-c(23))/c(24)/sqrt(2)).^2-intP*((x(:,2)-c(25))/c(26)/sqrt(2)).^2)+ 90
% c(27)*exp(-intP*((x(:,1)-c(28))/c(29)/sqrt(2)).^2-intP*((x(:,2)-c(30))/c(31)/sqrt(2)).^2)+ 90
% c(32)*exp(-intP*((x(:,1)-c(33))/c(34)/sqrt(2)).^2-intP*((x(:,2)-c(35))/c(36)/sqrt(2)).^2)+ 90
% c(37)*exp(-intP*((x(:,1)-c(38))/c(39)/sqrt(2)).^2-intP*((x(:,2)-c(40))/c(41)/sqrt(2)).^2)
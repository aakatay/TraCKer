load test_gausFit_highVar
i2=1;
%fun11var = @(c,x)c(1)+isAnotherSpot*c(2)*exp(-intPG*((x(:,1)-c(3))/c(4)/sqrt(2)).^2-intPG*((x(:,2)-c(5))/c(6)/sqrt(2)).^2)+isAnotherSpot*c(7)*exp(-intPG*((x(:,1)-c(8))/c(9)/sqrt(2)).^2-intPG*((x(:,2)-c(10))/c(11)/sqrt(2)).^2);
%funSpNonReg = fun11var
LB(11)=[];
UB(11)=[];
[fitVal_,errFit7by7(i),err2res_,EXITFLAG] = gausFit(data7by7diff,funSpNonReg,G0nonReg,LB(:,i2),UB(:,i2));

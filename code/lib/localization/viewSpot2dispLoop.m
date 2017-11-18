%%

%% exit flag
switch EXITFLAG
    case 1
        exitflag(i) = {'lsqcurvefit converged to a solution'};
    case 2
        exitflag(i) = {'Change in X too small'};
    case 3
        exitflag(i) = {'Change in RESNORM too small'};
    case 4
        exitflag(i) = {'Computed search direction too small'};
    case 0
        exitflag(i) = {'Too many function evaluations or iterations'};
    case -1
        exitflag(i) = {'Stopped by output/plot function'};
    case -2
        exitflag(i) = {'Bounds are inconsistent'};
end        

%  3by3 rectangle position
temp = zeros(7); temp(pxYc-1,pxXc-1)=1; temp = padarray(temp,[1 1]);
[yCoor3by3(i) xCoor3by3(i)] = find(temp==1); % rectangle coord for 3 by 3
     
% 2by2 peak
data2by2 = data3by3(linIx); 
data2by2 = reshape(data2by2,[2 2]);%gaussian reshaped as matrix

% 2by2 error
Ifit2_= funSpNonReg(fitVal,x); % gaussian fit result
Ifit2(:,:,i)=reshape(Ifit2_,[n m]);%gaussian reshaped as matrix
fit3by3 = Ifit2(pxYc-1:pxYc+1,pxXc-1:pxXc+1);
fit2by2 = reshape(fit3by3(linIx),2,2);
err2by2(i) = sqrt(mean(mean((data2by2-fit2by2).^2))).^(1/intP)*100;

% errors (pixel intensity equivalent )
errFit7by7(i) = (sqrt(errFit7by7(i))/numel(Window7by7) ).^(1/intP); %(=RESNORM) : sum {(FUN(X,XDATA)-YDATA).^2}
%     Irat(i)= min([fitVal(2,i) fitVal(3,i)])/max([fitVal(2,i) fitVal(7,i)]); % intensity ratios
%     dist(i) = sqrt((fitVal(3,i)-fitVal(8,i))^2 + (fitVal(5,i)-fitVal(10,i))^2); % distance betw. spots

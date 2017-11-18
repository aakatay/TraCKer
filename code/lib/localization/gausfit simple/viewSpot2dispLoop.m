%%
% tP = 0.5; % tolerance position            
c02G5by5([3,5]) = [px5Xmax,py5Ymax]; lb2G5by5([3,5]) = [px5Xmax,py5Ymax]-tP; ub2G5by5([3,5]) = [px5Xmax,py5Ymax]+tP;  
c02G5by5([8,10]) = [0,0]; lb2G5by5([8,10]) = [-inf,-inf]; ub2G5by5([8,10]) = [inf,inf];  
c02G5by5([1,2,7]) = [0 1 0];  lb2G5by5([1,2,7]) = [-minmaxIntTol 1/tol 1/tol];  ub2G5by5([1,2,7]) = [minmaxIntTol 1*tol 1*tol];  
is2ndDoubleFit = 0;
if is2ndDoubleFit
    [cc2(:,sp),err2G5by5(sp),err2res_,EXITFLAG(2)] = gausFit(Window5by5,fun2,c02G5by5,lb2G5by5,ub2G5by5);
else
    cc2(:,sp) = nan(11,1);
    err2G5by5(sp) = nan;
    EXITFLAG(2) = nan;
end
Ifit2_=fun2(cc2_2,x); %your fitted gaussian in vector
Ifit2(:,:,sp)=reshape(Ifit2_,[n m]);%gaussian reshaped as matrix

for j = 1: numel(EXITFLAG)
    switch EXITFLAG(j)
        case 1
            exitflag(sp,j) = {'lsqcurvefit converged to a solution'};
        case 2
            exitflag(sp,j) = {'Change in X too small'};
        case 3
            exitflag(sp,j) = {'Change in RESNORM too small'};
        case 4
            exitflag(sp,j) = {'Computed search direction too small'};
        case 0
            exitflag(sp,j) = {'Too many function evaluations or iterations'};
        case -1
            exitflag(sp,j) = {'Stopped by output/plot function'};
        case -2
            exitflag(sp,j) = {'Bounds are inconsistent'};
    end        
end            
 
err2G5by5_2(sp) = nan;
intMul = intMax(sp)*intPeak_(sp);
if isTryFit(sp)      
    % 2by2 peak
    %data2by2 = data3by3(linIx); 
    %data2by2 = reshape(data2by2,[2 2]);%gaussian reshaped as matrix
    %meanSpotPeak(sp) = mean(data2by2(:));     

     temp = zeros(5); temp(pxYc-1,pxXc-1)=1; temp = padarray(temp,[2 2]);
     [yCoor5by5(sp) xCoor5by5(sp)] = find(temp==1); % rectangle coord for 3 by 3
     
     fit3by3 = Ifit2(pxYc-1:pxYc+1,pxXc-1:pxXc+1,sp);           
%     fit2by2 = reshape(fit3by3(linIx),2,2);     
%     err2by2(sp) = sqrt(mean(mean((data2by2-fit2by2).^2))).^(1/intP)*100; 

     err3by3(sp) = sqrt(mean(mean((data3by3-fit3by3).^2))).^(1/intP)*100; 
     err2G5by5_2(sp) = err3by3(sp);
    
    Irat_2(sp)= min([cc2_2(2,sp) cc2_2(7,sp)])/max([cc2_2(2,sp) cc2_2(7,sp)]); % intensity ratios
    cc2disp_2(:,sp) = cc2_2(:,sp)+[0 0 2 0 2 0 0 2 0 2 0]';
    Ifit2_2_=fun2_2(cc2_2,x); %your fitted gaussian in vector
    Ifit2_2(:,:,sp)=reshape(Ifit2_2_,[n m]);%gaussian reshaped as matrix
    dist_2(sp) = sqrt((cc2_2(3,sp)-cc2_2(8,sp))^2 + (cc2_2(5,sp)-cc2_2(10,sp))^2); % distance betw. spots
    fit3by3_2 = Ifit2_2(pxYc-1:pxYc+1,pxXc-1:pxXc+1);
%     fit2by2_2 = reshape(fit3by3_2(linIx),2,2);
%     err2by2_2(sp) = sqrt(mean(mean((data2by2-fit2by2_2).^2))).^(1/intP2G)*100;
    Ifit2SHOW_2_=fun2_2(cc2disp_2(:,sp),xSHOW); %your fitted gaussian in vector
    Ifit2SHOW_2(:,:,sp)=reshape(Ifit2SHOW_2_,[n2 m2]);%gaussian reshaped as matrix
    Ifit2SHOW_2(:,:,sp) = Ifit2SHOW_2(:,:,sp)*intMul;
end

IdatSHOW(:,:,sp) = (spotWin(:,:,i)/intMax(sp))/intPeak;
err1G5by5(sp) = (sqrt(err1G5by5(sp))/numel(Window5by5) ).^(1/intP1G)*100; % RESNORM : sum {(FUN(X,XDATA)-YDATA).^2}
err2G5by5(sp) = (sqrt(err2G5by5(sp))/numel(Window5by5) ).^(1/intP)*100;
Irat(sp)= min([cc2(2,sp) cc2(7,sp)])/max([cc2(2,sp) cc2(7,sp)]); % intensity ratios

% cc1disp(:,sp) = cc1(:,sp)+[0 0 2 0 2 0]'; 
% cc2disp(:,sp) = cc2(:,sp)+[0 0 2 0 2 0 0 2 0 2 0]';
% Ifit1SHOW_=fun1(cc1disp(:,sp),xSHOW); %your fitted gaussian in vector
% Ifit1SHOW(:,:,sp)=reshape(Ifit1SHOW_,[n2 m2]);%gaussian reshaped as matrix            
% Ifit2SHOW_=fun2(cc2disp(:,sp),xSHOW); %your fitted gaussian in vector
% Ifit2SHOW(:,:,sp)=reshape(Ifit2SHOW_,[n2 m2]);%gaussian reshaped as matrix
% Ifit2SHOW(:,:,sp) = Ifit2SHOW(:,:,sp)*intMul;

IdatSHOW(:,:,sp) = IdatSHOW(:,:,sp)*intMul;


[py,px] = find(Idat(:,:,sp)==max(max(Idat(:,:,sp))));
[vv ix_] = sort((py-3).^2+(px-3).^2);
py = py(ix_(1));px = px(ix_(1));
temp = zeros(5); temp(py,px)=1; temp = padarray(temp,[2 2]);
[YpxPeak(sp) XpxPeak(sp)] = find(temp==1);
dist(sp) = sqrt((cc2(3,sp)-cc2(8,sp))^2 + (cc2(5,sp)-cc2(10,sp))^2); % distance betw. spots

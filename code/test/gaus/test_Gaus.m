% tests errors due to pixelation
clear
close all
xNorm = 40/100; yNorm = 40/100;
xNorm = 35/100;yNorm = 35/100;
xNorm = 50/100;yNorm = 50/100;

k = 1;
m = 5*k; n = 5*k; 
[X,Y]=meshgrid(1:n,1:m);%your x-y coordinates
x(:,1)=X(:); % x= first column
x(:,2)=Y(:); % y= second column
sg0 = 0.73; %[px] sg = 0.73nm (FWHM=183nm)


fun1 = @(c,x) c(1)+c(2)*exp(-((x(:,1)-c(3))/c(4)/sqrt(2)).^2-((x(:,2)-c(5))/c(6)/sqrt(2)).^2);
sg = sg0*k;
cx = m*xNorm+0.5;
cy = n*yNorm+0.5;
cc = double([0 1 cx sg cy sg]); %start-guess here

gaus = fun1(cc,x);
gaus = reshape(gaus,m,n);
gaus = gaus/max(gaus(:));

figure(1)
imagesc(gaus)


%%
clear x
k = 100;
m = 5*k; n = 5*k; 
[X,Y]=meshgrid(1:n,1:m);%your x-y coordinates
x(:,1)=X(:); % x= first column
x(:,2)=Y(:); % y= second 
fun1 = @(c,x) c(1)+c(2)*exp(-((x(:,1)-c(3))/c(4)/sqrt(2)).^2-((x(:,2)-c(5))/c(6)/sqrt(2)).^2);


%% change position here
xNorm(1) = 50/100;yNorm(1) = 50/100;
xNorm(2) = 55/100;yNorm(2) = 50/100;
xNorm(2) = 60/100;yNorm(2) = 50/100;
xNorm(2) = 60/100;yNorm(2) = 60/100;

%% gaussians
SNR = 2;
SNR = inf;
for r = 1:2
    sg = sg0*k;
    cx = m*xNorm(r)+0.5;
    cy = n*yNorm(r)+0.5;
    cc = double([0 1 cx sg cy sg]); %start-guess here

    gausHR_ = fun1(cc,x);
    gausHR(:,:,r) = reshape(gausHR_,m,n);
    gausHR(:,:,r) = gausHR(:,:,r)/max(max(gausHR(:,:,r)));

    %% integrate
    for i = 1:m/k
        for j = 1:n/k
            xx = (i-1)*k+1;
            yy = (j-1)*k+1;
            GAUS(j,i,r) = sum(sum(gausHR(yy:yy+k-1,xx:xx+k-1,r)));
        end
    end
end
    GAUS = GAUS/max(GAUS(:));
    GAUS = GAUS + 1/SNR;

figure(3)
imagesc(GAUS(:,:,1))

figure(4)
imagesc([gaus  ones(5,2) GAUS(:,:,1) ones(5,2) GAUS(:,:,1)-gaus]);
axis image
colorbar

sumInt1 = sum(sum(GAUS(2:4,2:4,1)));
sumInt2 = sum(sum(GAUS(2:4,2:4,2)));

peakInt1 = max(max(GAUS(2:4,2:4,1)));
peakInt2 = max(max(GAUS(2:4,2:4,2)));

sumInt2by2_1 = sum(sum(GAUS(3:4,3:4,1)));
sumInt2by2_2 = sum(sum(GAUS(3:4,3:4,2)));

figure(5)
subplot(2,1,1)
imagesc( [gausHR(:,:,1) ones(m,2*k) gausHR(:,:,2) ones(m,2*k) gausHR(:,:,1)-gausHR(:,:,2)]);
axis image
colorbar
subplot(2,1,2)
imagesc( [GAUS(:,:,1) ones(5,2) GAUS(:,:,2) ones(5,2) GAUS(:,:,1)-GAUS(:,:,2)]);
title([sprintf('intensity1:%.02f, intensity2by2:%.02f, peak1:%.02f',sumInt1,sumInt2by2_1,peakInt1);...
    sprintf('intensity2:%.02f, intensity2by2:%.02f, peak2:%.02f',sumInt2,sumInt2by2_2,peakInt2)])
axis image
colorbar


%% results
display('in case the peak intensity is spread to four pixels total intensity in 3by3 is 9.75% lower, peak int. is 1/3 lower')
display('due to possible neighbouring spots in tw0 pixels distance only only peak intensity should be used as the starting guess')
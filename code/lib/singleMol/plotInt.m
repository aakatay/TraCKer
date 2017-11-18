if ~exist('intSum.mat')
    fname = 'acq031_-X-Y376x154_279X8Y12x12_1-26004_crop1.tif';
    imginfo = imfinfo(fname);
    N = length(imginfo);
    for i = 1:N
        A(:,:,i) = imread(fname,i);
    end
    intSum = squeeze(mean(mean(A,1),2));
    save('intSum','intSum')
else
    load('intSum')
end


%%
n1= 1011;
intSum = intSum(n1:end);

intSumBckGrnd = walkingMin(intSum,500);
figure(5)
plot(intSum); hold; plot(intSumBckGrnd,'r');

N = length(intSumBckGrnd);
x=1:N;
plot(intSum);
a=intSum(1);
b=-0.0001;
c=intSum(1);
d=0;    

bleach1 = fit(x',intSumBckGrnd','exp2','StartPoint',[a b c d]);

a=bleach1.a;
b=bleach1.b;
c=bleach1.c;
d=bleach1.d;

% calc the fit background
N = length(intSum);
x=1:N;
bleachFit = a*exp(b*x) + c*exp(d*x);
bleachFit = bleachFit';

figure(1)
plot(bleach1,x,intSum)
grid minor

figure(2)
plot(bleachFit)


intSumCorr = intSum-bleachFit - min(intSum-bleachFit);

figure(3)
plot([1:N]+n1 ,intSumCorr);
grid minor


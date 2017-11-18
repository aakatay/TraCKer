clear all;
close all;
load('intSum.mat')

intSum = intSum(1011:end);

N = length(intSum);
x=1:N;
plot(intSum);
a=intSum(1);
b=-0.0001;
exp1 = a*exp(b*x);
hold;plot(exp1,'r');hold;
grid minor

plot(exp1,x,intSum')
bleach1 = fit(x',intSum,'exp2','StartPoint',[intSum(1) -0.001 intSum(1) 0])

a=bleach1.a;
b=bleach1.b;
c=bleach1.c;
d=bleach1.d;

bleachFit = a*exp(b*x) + c*exp(d*x);

plot(bleach1,x,intSum)

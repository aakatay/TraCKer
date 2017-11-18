clear all;
close all;
x = 1:5;
[X Y] = meshgrid(x,x);

cx = 3.5;
cy = 1.5;
sg = 10;
I = exp(-((X-cx)/sg/sqrt(2)).^2-((Y-cy)/sg/sqrt(2)).^2);
fun = @(c,x) c(1)+c(2)*exp(-((x(:,1)-c(3))/c(4)/sqrt(2)).^2-((x(:,2)-c(5))/c(4)/sqrt(2)).^2);

c0 = [0 1 2.5 1.5 2.5];
cc = gausFit(I,fun,c0)
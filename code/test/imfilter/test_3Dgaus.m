% not used

N = 3;
x = -N:N;
X = ndgrid(x,x,x);
%X1 = permute(X,[]);
X2 = permute(X,[3 2 1]);
X3 = permute(X,[2 1 3]);

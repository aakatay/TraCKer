% sparse

SZ = 1e5; % matrix size
N = 1e3; % number of elements

x = floor(rand(N,1)*SZ);
y = floor(rand(N,1)*SZ);
z = uint16(x); 


A = sparse()
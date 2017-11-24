x = ceil(rand(1000,1)*512);
y = ceil(rand(1000,1)*512);
szYX=[512 512];
szYXF=[512 512 1000];

A = zeros(512);
tic; A(sub2ind(szYX,x,y)) = 1; toc;
A = zeros(512,512,1000);
tic; A(sub2ind(szYXF,x,y,[1:1000]')) = 1; toc;
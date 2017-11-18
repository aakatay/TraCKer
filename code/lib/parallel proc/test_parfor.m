clear
N = 1e3;
sz = 1e3;
img = rand(sz);

tic 
for i = 1:N 
    img = rand(sz);
    A(:,:,i) = img.^2;
end
toc

tic
parfor i=1:N
    img = rand(sz);
    A(:,:,i) = img.^2;
end
toc


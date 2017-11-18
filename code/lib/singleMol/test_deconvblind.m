

a0 = rand(100,100)*1.01;
a = double(im2bw(a0,1)); 
% PSF
gausKernelSz = 20;
gausKernelSg = 4;
PSF = fspecial('gaussian', gausKernelSz, gausKernelSg);

A = conv2(a,PSF,'same');
DAMPAR=A/100;
[a2,~] = deconvblind(A,PSF,20,DAMPAR);

subplot(2,2,1)
imagesc(a)
subplot(2,2,2)
imagesc(A)
subplot(2,2,3)
imagesc(a2)
subplot(2,2,4)
imagesc(a-a2)







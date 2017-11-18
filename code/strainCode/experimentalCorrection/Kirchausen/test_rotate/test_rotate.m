clear
img = zeros(100);
img(71:100,31:90)=1;
rotAngle = -10;
imgRot = imrotate(img,rotAngle,'bicubic');
figure(1)
n=4;
subplot(1,n,1)
imagesc(img)
axis equal
axis tight
subplot(1,n,2)
imagesc(imgRot)
axis equal
axis tight
xlabel('z')
ylabel('x')
title('deskewed (isotropic sampling)')

imgRotLowSamp = imgRot(:,1:4:size(imgRot,2));
subplot(1,n,3)
imagesc(imgRotLowSamp)
xlabel('z''=4*z')
ylabel('x')
title('deskewed (nonisotropic sampling ie. z''=4*z)')
axis equal
axis tight

% shear correction (non isotropic rotation)

rotAngle2=29;
subplot(1,n,4)
imgRotLowSampBackRotated = imrotate(imgRotLowSamp,rotAngle2,'bicubic')
imagesc(imgRotLowSampBackRotated)
axis equal
axis tight
maximize


return;

sy= 0;%1/tan(rotAngle/180*pi); % tan theta;
sx = tan(rotAngle/180*pi); % tan theta
xx = 1;
yy = 1;
zz = 1;
xform =  [xx sx 0; sy yy 0; 0 0 zz];
tform = affine2d(xform);
imgShear1 = imwarp(imgRotLowSamp,tform);
imagesc(imgShear1)

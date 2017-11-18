
clear
img = zeros(100,214);
img(71:100,6:165)=repmat([1:160],30,1);
img(65:70,6:165)=160;
img(61:65,6:165)=0;
img(55:60,6:165)=160;
figure(1)
maximize
n=5;

%% 1
subplot(1,n,1)
imagesc(img)
axis equal
axis tight
title('real')
%% shearing
shear = 1.0214;
sx = shear;
sy = 0; % tan theta
xx = 1;
yy = 1;
xform =  [xx sy 0; sx yy 0; 0 0 1];
tform = affine2d(xform)
imgShear_ = imwarp(img,tform,'cubic');
imgShear = zeros(380,size(imgShear_,2));
imgShear(1:size(imgShear_,1),:) = imgShear_
%% 2
subplot(1,n,2)
imagesc(imgShear)
axis equal
axis tight
xlabel('z')
ylabel('x')
title('sheared (isotropic sampling)')
%% noniso sampling shear image
imgShearLowSamp = imgShear(:,1:4:size(imgShear,2));
%% 3
subplot(1,n,3)
imagesc(imgShearLowSamp)
axis equal
axis tight
xlabel('z')
ylabel('x')
title('sheared (nonisotropic sampling ie. z''=4*z)')
%% deshearing
sx = -shear;
sy = 0; % tan theta
xx = 4;
yy = 1;
xform =  [xx sy 0; sx yy 0; 0 0 1];
tform = affine2d(xform)
imgDeShear = imwarp(imgShearLowSamp,tform);

%% 4
subplot(1,n,4)
imagesc(imgDeShear)
axis equal
axis tight
xlabel('z')
ylabel('x')
title('desheared (nonisotropic sampling ie. z''=4*z)')
% calculating the theta for non iso
x=4; y=1;
thISO = atan(shear)/pi*180;
thNONISO = atan(tan(thISO/180*pi)*x/y)/pi*180
% deshearing low sampled img
imgShearLowSamp

%% 5
subplot(1,n,5)
imagesc(imgDeShear)
axis equal
axis tight
xlabel('z')
ylabel('x')
title('desheared (nonisotropic sampling ie. z''=4*z)')



% find the shearing angle :thISO 
% th : shearing of the acquired image (noniso sampled)
px = 0.104;
pz = 0.38;
th = 18.5;

thISO = atan(tan(th/180*pi)*px/pz)/pi*180

return;
% generate skew data
clear
a = zeros(30,100);

sVec = 0:5:495;
for i = 1:size(a,2)
    a(i,:) = sVec+i;
end

sx = 0;
sy = .2; % tan theta
xform =  [1 sx 0; sy 1 0; 0 0 1];
tform = affine2d(xform)
J = imwarp(a,tform);
figure
subplot(1,2,1)
imagesc(a)
subplot(1,2,2)
imagesc(J)

imwrite(uint16(a),'shearedImg.tif');
%stackWrite(a,'shearedStack.tif');


















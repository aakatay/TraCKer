
cd('E:\MATLAB\image type')


a=logical(round(rand(128,116,100)*2));

%a = uint16(a);
stackWrite(a,'logical.tif',1)
% error : cannot open tiff files compressed in this fashion (2)



%% 10 by 10 - 2 frame

a = logical(round(rand(10,10,2)*2));
stackWrite(a,'logical_10by10.tif')
% error : cannot open tiff files compressed in this fashion (2)

%% 10 by 10 - 1 frame

a = logical(round(rand(10,10,1)*2));
stackWrite(a,'logical_10by10.tif')
% error : cannot open tiff files compressed in this fashion (2)


%% 10 by 10 - 2 frame - 16bit

a = logical(round(rand(10,10,2)*2));
a = uint16(a);
stackWrite(a,'logical_10by10.tif')

% works as 16 bit

%% 10 by 10 - 1 frame

a = logical(round(rand(10,10,1)*2));
imwrite(a,'logical_10by10.tif')
% error : cannot open tiff files compressed in this fashion (2)

%% 10 by 10 - 1 frame - compression

a = logical(round(rand(10,10,1)*2));
imwrite(a,'logical_10by10.tif', 'Compression','packbits')

% works as (1 bit), but 1 frame only

%% 10 by 10 - 2 frame - compression

a = logical(round(rand(10,10,2)*2));
imwrite(a(:,:,1),'logical_10by10.tif', 'Compression','packbits')
imwrite(a(:,:,2),'logical_10by10.tif', 'Compression','packbits','WriteMode','append')

% displays only 1 frame (1 bit), should be 2 frames

%% 10 by 10 - 2 frame - compression

a = logical(round(rand(10,10,2)*2));
imwrite(a(:,:,1),'logical_10by10.tif', 'Compression','packbits')
imwrite(a(:,:,2),'logical_10by10.tif', 'Compression','packbits','WriteMode','append')

% displays only 1 frame (1 bit), should be 2 frames

%%

imfinfo('logical_10by10.tif')
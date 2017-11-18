cd C:\MATLAB\TraCKer\data_eclipse\EMCCD
clear
close

d = rdir('*.tif');
[Y,J] = sort([d.datenum]); 

for i = 1:size(d,1)
    fn(i) = {d(J(i)).name}; % filename
end
imagesc(imread(fn{10}));

if ~exist('XY.mat')
    h = figure(1);
    btn = 0;
    [x y] = ginput;
    save('XY','x','y');
else
    load XY
end


nS = size(x,1);
D = 8; % spot window half  width
for i = 1:nS
    x1(i)= x(i)-D; x2(i)= x(i)+D;
    y1(i)= y(i)-D; y2(i)= y(i)+D;
end
x1 = round(x1);
x2 = round(x2);
y1 = round(y1);
y2 = round(y2);

nFrm = 95; % number of frames
nAcq = 19; % number of acquistion sets
%S = zeros()
for s = 1:nS
    for i = 1:nAcq
        for j = 1:nFrm
            A(:,:,(i-1)*nFrm+j) = imread(cell2mat(fn(i)),j);    
        end
    end
    S(:,:,:,s) = A(y1(s):y2(s),x1(s):x2(s),:);
end
save('Spots','S');


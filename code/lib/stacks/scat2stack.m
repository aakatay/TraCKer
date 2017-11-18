%converts positional info to image info 
clear;
dataXY = dir('posData*.mat');
load(dataXY.name);
frames = size(X,2);

data = dir('data*.tif');
a = imread(data.name,1);


% determine magnification
mag = 120;
sizeNotOK = 1;    
while sizeNotOK    
    v=size(a).*mag;
    m=v(2);
    n=v(1);
    pos=get(0,'ScreenSize');
    pos=pos(3:4) - [m n];
    if isempty(find(pos<0))
        sizeNotOK = 0;
    else
        mag= mag-1;
    end
end    

X = ceil(mag*X);
Y = ceil(mag*Y);

ix = [];
for i = 1:frames
    if isempty(find(X(:,i)~=0)), continue; end;
    ix_ = sub2ind([m,n,frames],X(find(X(:,i)~=0),i),Y(find(Y(:,i)~=0),i),repmat(i,numel(find(Y(:,i)~=0)),1));
    ix = [ix; ix_];
end

imgScat = zeros(n*m*frames,1);
imgScat(ix) = 1;
imgScat = reshape(imgScat,[n,m,frames]);



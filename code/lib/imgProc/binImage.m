function img = binImage(img,binsize)
    m = binsize;
    n1=size(img,1)/m; n2=size(img,2)/m;
    if n1~=round(n1) || n2~=round(n2), error('size error'); end;
    c = reshape(img,[m n1 m n2]);
    c=sum(c,1);
    c=sum(c,3);
    img=reshape(c,[n1 n2]);
end
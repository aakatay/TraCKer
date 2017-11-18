clear all
f=rdir('*.tif');
N = size(f,1);
for i = 1:N
    nm = f(i).name;
    %tnm(i) = str2num(nm(47:53));
    tnm(i) = str2num(nm(27:33));
    
end
TNM = tnm;
diff=TNM(2:end)-TNM(1:end-1);
diff'

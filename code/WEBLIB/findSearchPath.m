clear all;

mP = matlabpath;
nP = numel(mP);

srchStr='Coeff.mat';
iLast = 0;
for i = 1: nP
    if mP(i) == ';'
        filePath{1} = mP(iLast+1:i-1);
        iLast = i;        
    end
end

for i= 1:numel(filePath)
    a = rdir(cell2mat(filePath(i)),'strfind(name, ''Coeff.mat'')');
    if numel(a)
        disp(filePath(i));  
    end
end

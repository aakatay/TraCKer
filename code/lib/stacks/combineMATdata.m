clear
fname = dir('posData*');

XX = []; YY = [];
for i = 1:numel(fname)
    fn = fname(i).name;
    load(fn)
    XX = [XX; X(:)];
    YY = [YY; Y(:)];
end
X = XX;
Y = YY;
save('posData-Combine','X','Y')

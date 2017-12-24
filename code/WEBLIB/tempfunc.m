clear
fname = dir('posData*');

for i = 1:numel(fname)
    fn = fname(i).name;
    load(fn)
    XX = X;
    YY = Y;
end

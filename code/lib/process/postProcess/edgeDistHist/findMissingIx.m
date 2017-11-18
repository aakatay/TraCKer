function ix = findMissingIx(ix1in,ix2in)
if numel(ix1in)>numel(ix2in)
    ix1 = ix1in; ix2 = ix2in; ixtxt='ix2';
else
    ix1 = ix2in; ix2 = ix1in; ixtxt='ix1';
end

for i = 1:numel(ix1)
    if isempty(find(ix2==ix1(i)))
        disp(sprintf('missing ix: _%02i in %s',ix1(i),ixtxt));
    end
end



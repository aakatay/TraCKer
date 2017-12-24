function [mx mn] = dispRange(A)
    mx = max(A(:));
    mn = min(A(:));
    display(sprintf('max:%f, min:%f',mx,mn));
end
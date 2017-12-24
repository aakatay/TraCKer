function A = invert(A);
% invert BW color

mx = max(A(:));
A = mx-A+1;

end
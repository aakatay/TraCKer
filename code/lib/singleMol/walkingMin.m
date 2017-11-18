function B = walkingMin(A,winSz)
% returns the walking minimum of the vector
    N = length(A);
    for i = 1:N-winSz+1
        B(i) = min(A(i:i+winSz-1));
    end
        
end
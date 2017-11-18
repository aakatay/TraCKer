% finds the positions of the nonzero elements in the array 
% INT: # of occurances
function [Y,X,INT] = findPos(A)
    INT = A(A>0);
    [Y,X]=find(A>0);
end

%%
% finds the positions of the nonzero elements in the array such that:
% number of localizations = total intensity    INT = A(A>0);
%     [Y,X]=find(A>0);
%     N =  numel(INT);
%     for i = 1:N
%         int = INT(i)-1;
%         Y = [Y; repmat(Y(i),int,1)];
%         X = [X; repmat(X(i),int,1)];
%     end
%     cc=3;
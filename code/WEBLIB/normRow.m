function [ y] = normRow(x)
%based on NORMR Normalize rows of matrices.

% Checks
if nargin < 1,nnerr.throw('Not enough input arguments.'); end
wasMatrix = ~iscell(x);
x = nntype.data('format',x,'Data');

% Compute
y = cell(size(x));
for i=1:numel(x)
  xi = x{i};
  cols = size(xi,2);
  n = 1 ./ sqrt(sum(xi.*xi,2));
  yi = xi .* n(:,ones(1,cols));
  yi(~isfinite(yi)) = 0;
  y{i} = yi;
end

% Format
if wasMatrix, y = y{1}; end

    
end 
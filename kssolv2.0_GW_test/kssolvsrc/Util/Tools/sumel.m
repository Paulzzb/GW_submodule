function s = sumel(A)
% SUMEL sums all the elements.
%
%   s = SUMEL(A) returns the total sum of all elementals in A.
%
%   See also sum.

if iscell(A)
    s = sumel(cellfun(@sumel,A));
else
    s = sum(A(:));
end

end

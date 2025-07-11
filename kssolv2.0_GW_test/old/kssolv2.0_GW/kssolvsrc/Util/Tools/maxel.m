function m = maxel(A)
% MAXEL maximum of all the elements.
%
%   m = MAXEL(A) returns the maximum of all elementals in A.
%
%   See also sumel, minel.

if iscell(A)
    m = maxel(cellfun(@maxel,A));
else
    m = max(A(:));
end

end

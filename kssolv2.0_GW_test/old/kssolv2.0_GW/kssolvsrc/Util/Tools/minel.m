function m = minel(A)
% MINEL minimum of all the elements.
%
%   m = MINEL(A) returns the minimum of all elementals in A.
%
%   See also sumel, maxel.

if iscell(A)
    m = minel(cellfun(@minel,A));
else
    m = min(A(:));
end

end

function b = iszero(x)
% ISZERO Check value x.
%
%   b = ISZERO(x) returns whether x is empty or zero.
%
%   See also isempty.

b = isempty(x) || (sumel(x ~= 0) == 0);

end

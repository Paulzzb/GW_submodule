function k = doublefactorial(n)
% DOUBLEFACTORIAL double facotorial function.
%   k = DOUBLEFACTORIAL(n) for integer n, is the double factorial of n,
%   i.e., k = n*(n-2)*(n-4)*...*2 for even n, and k = n*(n-2)*...*1 for odd
%   n.
%
%   See also factorial.

k = prod(n:-2:1);

end

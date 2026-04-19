% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08

function d = divisors_int32(n)
% Positive divisors of scalar integer n (column int32).

  n = int32(n);
  if n < 1
    d = int32([]);
    return;
  end

  d = int32([]);
  for k = int32(1):n
    if mod(n, k) == 0
      d(end + 1, 1) = k; %#ok<AGROW>
    end
  end
end

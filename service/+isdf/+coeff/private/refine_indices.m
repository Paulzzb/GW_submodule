% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/07 ZZ
%
% Ported from kssolvGW .../ISDF/refine_indices.m

function indices = refine_indices(indices, weight, rpoint, tol)
%REFINE_INDICES  Replace invalid / tiny-weight centroid indices.

  Npoint = length(weight);
  Ncenter = length(indices);

  for i = 1:Ncenter
    current = indices(i);
    if current < tol
      ref = rpoint(i, :);
      diff = rpoint - ref;
      dist = sqrt(sum(diff.^2, 2));
      [~, out] = sort(dist, 'ascend');
      for ip = 2:Npoint - 1
        if weight(out(ip)) > tol
          indices(i) = out(ip);
          break
        end
      end
    end
  end
end

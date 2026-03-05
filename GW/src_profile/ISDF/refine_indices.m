function indices = refine_indices(indices, weight, rpoint, tol)
  % Replace indices with indices that are within a tolerance of the
  %
  Npoint = length(weight);
  Ncenter = length(indices);
  %
  for i = 1:Ncenter
    current = indices(i);
    if current < tol
      % Find a point near that point to replace it
      ref = rpoint(i, :);
      diff = rpoint - ref;
      dist = sqrt(sum(diff.^2, 2));
      [~, out] = sort(dist, 'ascend');
      for ip = 2:Npoint-1
        if (weight(ip) > tol)
          indices(i) = out(ip);
          break
        end
      end
    end
  end

  
end % function
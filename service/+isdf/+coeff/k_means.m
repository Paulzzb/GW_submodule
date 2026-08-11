% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/07 ZZ
%
% Ported from kssolvGW .../ISDF/k_means.m

function ind_mu = k_means(rk, weight, opt)
%K_MEANS  Weighted k-means on the fine FFT grid (Cartesian points).
%
%   ind_mu = isdf.coeff.k_means(rk, weight, opt)
%
%   opt.seed       - RNG seed
%   opt.init       - 'sequence'|'random'|'wrs'|'kmeans++'|'kmeans++_false'|'atom'
%   opt.fftgrid    - [n1 n2 n3]
%   opt.supercell  - 3x3 (same layout as legacy sys.supercell / d_lat.a1a2a3)
%   opt.atompos    - nat x 3 Cartesian (for init='atom'); may be empty

  weight_tol = 1e-5;

  n1 = double(opt.fftgrid(1));
  n2 = double(opt.fftgrid(2));
  n3 = double(opt.fftgrid(3));
  [I, J, K] = ndgrid( ...
    (0:n1-1) / n1 - ((0:n1-1) >= n1 / 2), ...
    (0:n2-1) / n2 - ((0:n2-1) >= n2 / 2), ...
    (0:n3-1) / n3 - ((0:n3-1) >= n3 / 2));
  points = reshape(cat(4, I, J, K), [], 3) * double(opt.supercell);

  atompos = [];
  if isfield(opt, 'atompos') && ~isempty(opt.atompos)
    atompos = double(opt.atompos);
    if ndims(atompos) == 3
      atompos = reshape(atompos, size(atompos, 1), 3);
    end
  end

  [Nr, ~] = size(points);
  if rk > Nr
    output.err('The number of points (%d) is less than rk (%d).', Nr, rk);
  end

  seed = 0;
  if isfield(opt, 'seed') && ~isempty(opt.seed)
    seed = double(opt.seed);
  end
  rng(seed, 'twister');

  ind_mu = zeros(rk, 1);
  dist = zeros(Nr, rk);
  init = 'random';
  if isfield(opt, 'init') && ~isempty(opt.init)
    init = lower(char(string(opt.init)));
  end

  switch init
    case 'sequence'
      lastCentroids_ind_mu = round(1:floor(Nr / rk):Nr);
      lastCentroids_ind_mu = lastCentroids_ind_mu(1:rk)';

    case 'random'
      lastCentroids_ind_mu = randperm(Nr, rk)';

    case 'wrs'
      weight = abs(weight);
      lastCentroids_ind_mu = datasample(1:Nr, rk, 'Replace', false, 'Weights', weight);

    case 'kmeans++_false'
      lastCentroids_ind_mu = zeros(rk, 1);
      lastCentroids_ind_mu(1) = randi([1, Nr], 1);
      for i = 2:rk
        dist(:, i-1) = sum((points - points(lastCentroids_ind_mu(i-1), :)).^2, 2);
        min_dist = min(dist(:, 1:i-1), [], 2);
        weighted_dist = min_dist .* weight;
        lastCentroids_ind_mu(i) = randsample(1:Nr, 1, true, weighted_dist);
      end

    case 'kmeans++'
      lastCentroids_ind_mu = zeros(rk, 1);
      lastCentroids_ind_mu(1) = randi([1, Nr], 1);
      for i = 2:rk
        dist(:, i-1) = sum((points - points(lastCentroids_ind_mu(i-1), :)).^2, 2);
        min_dist = min(dist(:, 1:i-1), [], 2);
        weighted_dist = min_dist .* weight;
        [sorted_dist, Iord] = sort(weighted_dist, 'descend');
        rand_dist = 0.9 * rand * sum(sorted_dist);
        sum_dist = 0;
        for j = 1:Nr - 1
          sum_dist = sum_dist + sorted_dist(j);
          if sum_dist >= rand_dist
            lastCentroids_ind_mu(i) = Iord(j);
            break
          end
        end
      end

    case 'atom'
      if isempty(atompos)
        output.err('k_means init=''atom'' requires nonempty atom positions.');
      end
      lastCentroids_ind_mu = zeros(rk, 1);
      nat = size(atompos, 1);
      for i = 1:min(nat, rk)
        dist(:, i) = sum((points - atompos(i, :)).^2, 2);
        [~, Iord] = min(dist(:, i));
        lastCentroids_ind_mu(i) = Iord;
      end
      for i = nat + 1:rk
        dist(:, i-1) = sum((points - points(lastCentroids_ind_mu(i-1), :)).^2, 2);
        min_dist = min(dist(:, 1:i-1), [], 2);
        weighted_dist = min_dist .* weight;
        [sorted_dist, Iord] = sort(weighted_dist, 'descend');
        rand_dist = 0.9 * rand * sum(sorted_dist);
        sum_dist = 0;
        for j = 1:Nr - 1
          sum_dist = sum_dist + sorted_dist(j);
          if sum_dist >= rand_dist
            lastCentroids_ind_mu(i) = Iord(j);
            break
          end
        end
      end

    otherwise
      output.err('Unknown k-means initialization method: %s', init);
  end

  newCentroids_ind_mu = zeros(rk, 1);
  max_iteration = 100;
  iteration = 1;

  while true
    for nk = 1:rk
      dist(:, nk) = sum((points - points(lastCentroids_ind_mu(nk), :)).^2, 2);
    end
    [~, index_min] = min(dist, [], 2);
    shifted_points = points;

    center_new = newCentroids_ind_mu;
    for mk = 1:rk
      cluster = find(index_min == mk);
      if isempty(cluster)
        continue;
      end
      total_pointsMutilweight = sum(shifted_points(cluster, :) .* weight(cluster), 1);
      total_weight = sum(weight(cluster), 1);
      Centroids = total_pointsMutilweight / total_weight;
      distPoint2Centroid = sum((Centroids - shifted_points(cluster, :)).^2, 2);
      [~, Iord] = min(distPoint2Centroid);
      center_new(mk, 1) = cluster(Iord);
    end
    for mk = 1:rk
      center = newCentroids_ind_mu(mk);
      cluster = find(index_min == mk);
      if isempty(cluster)
        if center == 0
          t = 1;
        else
          t = center;
        end
        center_new(mk) = t;
      end
    end
    newCentroids_ind_mu = center_new;

    if max_iteration ~= 0 && iteration == max_iteration
      output.warn('K-Means reached max iterations (%d).', max_iteration);
      ind_mu = newCentroids_ind_mu;
      break
    elseif norm(newCentroids_ind_mu - lastCentroids_ind_mu) < 0.1
      ind_mu = newCentroids_ind_mu;
      break
    end

    iteration = iteration + 1;
    lastCentroids_ind_mu = newCentroids_ind_mu;
  end

  ind_mu = refine_indices(ind_mu, weight, points, weight_tol);
end

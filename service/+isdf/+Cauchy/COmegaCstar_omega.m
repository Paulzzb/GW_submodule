% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/06 ZZ

function [results, time, relError, iter] = COmegaCstar_omega(Phi, Psi, evOcc, evUnocc, options, omega)
% COMEGACSTAR  PhiPsi * Omega^{-1} * PhiPsi' via Cauchy integral.
%
%   Omega_ij = evOcc_i - evUnocc_j
%
%   Phi: n x n1 (occupied), Psi: n x n2 (unoccupied)
%   options.froErr, options.MaxIter

  if nargin < 6
    omega = 0;
  else
    if abs(real(omega)) > 1e-13
      output.err('real(omega) ~= 0.');
    end
  end  

  froErr = options.froErr;
  MaxIter = options.MaxIter;
  Error_array = [0, 0];
  ii = sqrt(-1);
  n = size(Phi, 1);
  M = evUnocc(end) - evOcc(end);
  m = evUnocc(1) - evOcc(end);
  if m <= 0
    output.err('No band gap in current system, Cauchy integral not applicable.');
  end
  k = (sqrt(M / m) - 1) / (sqrt(M / m) + 1);
  % [K, ~] = isdf.Cauchy.ellipk(k);
  % K(k): Driscoll ellipkkp(L) with L = -log(k)/pi  (same L as ellipjc below)
  L0 = -log(k) / pi;
  [K, ~] = isdf.Cauchy.ellipkkp(L0);


  results = zeros(n, n);
  discretePoints = [-K, K];
  newDiscretePoints = discretePoints;
  Iker = @(s) 1 ./ sqrt((1 + s.^2) .* (1 + k^2 .* s.^2));
  I = 1 / 2 * (integral(Iker, 0, k^(-1)));
  newGap = [];
  results_ = [];
  Error = inf;

  startTime = tic;

  for iter = 1:MaxIter
    oldResults = results;
    results = results / 2;

    for t = newDiscretePoints + I * ii
      [lambda, dlambda] = local_integrand(t, k, m, M);
      lambda = lambda + evOcc(end);
      OmegaOcc = diag(1 ./ (lambda - evOcc + omega));
      OmegaUnocc = diag(1 ./ (lambda - evUnocc));
      OccMatrix = Phi * OmegaOcc * Phi';
      UnoccMatrix = conj(Psi) * OmegaUnocc * conj(Psi)';
      matrix = (OccMatrix .* UnoccMatrix) .* (dlambda / pi / ii);

      if iter == 1
        results = results + K * real(matrix);
      else
        results = results + newGap * real(matrix);
      end
    end

    if iter >= 3
      oldResults_ = results_;
    end
    if iter >= 2
      results_ = (4 * results - oldResults) ./ 3;
    end

    if iter >= 3
      Error = norm(oldResults_ - results_, 'fro') / norm(results_, 'fro');
      Error_array(iter) = Error; 
      if Error <= froErr
        time = toc(startTime);
        results = results_;
        relError = Error;
        return
      end
    end

    newGap = (discretePoints(2) - discretePoints(1)) / 2;
    newDiscretePoints = discretePoints(1:end-1) + newGap;
    oldDiscretePoints = discretePoints;
    discretePoints = zeros(1, 2^iter + 1);
    discretePoints(1:2:end) = oldDiscretePoints;
    discretePoints(2:2:end) = newDiscretePoints;
  end

  output.warn( ...
    'Cauchy integral not converged after %d iterations (relerr=%.3e, tol=%.3e).', ...
    MaxIter, Error, froErr);
  relError = Error;
  results = results_;
  time = toc(startTime);
end

function [lambda, dlambda] = local_integrand(t, k, m, M)
  L = -log(k) / pi;
  [SN, CN, DN] = isdf.Cauchy.ellipjc(t, L);
  lambda = sqrt(m * M) .* ((k^(-1) + SN) ./ (k^(-1) - SN));
  dlambda = CN .* DN .* sqrt(m * M) .* ((2 * k.^(-1)) ./ (k^(-1) - SN).^2);
end

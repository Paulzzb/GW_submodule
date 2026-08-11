% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/29 ZZ

function [Q, R, info] = iterated_cholqr(X, tol, maxit)
%ITERATED_CHOLQR  Iterated Cholesky QR: X = Q*R with shift when chol fails.
%
%   [Q, R] = isdf.iterated_cholqr(X)
%   [Q, R, info] = isdf.iterated_cholqr(X, tol, maxit)
%
%   Implements Algorithm 4.1 (iterated CholeskyQR with shifts):
%     1. Q := X, R := I
%     2. A := Q'Q; R_tilde := chol(A) or chol(A + s*I) if chol breaks down
%        s = 11*(m*n + n*(n+1)) * u * ||X||_2^2
%     3. Q := Q / R_tilde,  R := R_tilde * R
%     4. Repeat until ||Q'Q - I||_F <= sqrt(n)*u
%
%   Inputs
%     X     鈥?m-by-n matrix (typically m >= n)
%     tol   鈥?optional stopping tolerance on ||Q'Q - I||_F (default sqrt(n)*eps)
%     maxit 鈥?optional maximum iterations (default 10)
%
%   Outputs
%     Q     鈥?m-by-n with approximately orthonormal columns
%     R     鈥?n-by-n upper triangular, X 鈮?Q*R
%     info  鈥?struct: .iters, .fro_err, .shift_count, .x_norm2

  if nargin < 1 || isempty(X)
    error('isdf:iterated_cholqr:Input', 'X must be a nonempty numeric matrix.');
  end
  if ~isnumeric(X) || ~ismatrix(X)
    error('isdf:iterated_cholqr:Input', 'X must be a 2-D numeric matrix.');
  end

  [m, n] = size(X);
  if n == 0
    Q = zeros(m, 0, class(X));
    R = zeros(0, 0, class(X));
    if nargout >= 3
      info = struct('iters', 0, 'fro_err', 0, 'shift_count', 0, 'x_norm2', 0);
    end
    return
  end
  if m < n
    error('isdf:iterated_cholqr:Size', ...
      'X must have at least as many rows as columns (m >= n); got %d-by-%d.', m, n);
  end

  u = eps(class(X));
  if nargin < 2 || isempty(tol)
    tol = sqrt(n) * u;
  else
    tol = double(tol);
  end
  if nargin < 3 || isempty(maxit)
    maxit = 10;
  else
    maxit = max(1, round(double(maxit)));
  end

  X = double(X);
  x_norm2 = norm(X, 2);
  if ~isfinite(x_norm2) || x_norm2 == 0
    x_norm2 = norm(X, 'fro');
  end

  Q = X;
  R = eye(n);
  shift_count = 0;
  fro_err = inf;

  for it = 1:maxit
    A = Q' * Q;
    A = (A + A') / 2;

    [R_tilde, used_shift] = local_chol_upper(A, m, n, u, x_norm2);
    shift_count = shift_count + double(used_shift);

    Q = Q / R_tilde;
    R = R_tilde * R;

    fro_err = norm(Q' * Q - eye(n), 'fro');
    if fro_err <= tol
      break
    end
  end

  if nargout >= 3
    info = struct( ...
      'iters', it, ...
      'fro_err', fro_err, ...
      'shift_count', shift_count, ...
      'x_norm2', x_norm2, ...
      'tol', tol);
  end
end

function [R_upper, used_shift] = local_chol_upper(A, m, n, u, x_norm2)
  used_shift = false;
  try
    R_upper = chol(A, 'upper');
  catch
    s = 11 * (m * n + n * (n + 1)) * u * (x_norm2^2);
    if ~isfinite(s) || s <= 0
      s = 11 * (m * n + n * (n + 1)) * u * max(norm(A, 'fro')^2, eps);
    end
    R_upper = chol(A + s * eye(n), 'upper');
    used_shift = true;
  end
end

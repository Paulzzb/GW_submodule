function [R, QQ_T, info] = safe_chol(A, m, tol, maxit)
%SAFE_CHOL  Iterated Cholesky QR on Gram matrix A = X'*X (Algorithm 4.1).
%
%   [R, QQ_T] = isdf.safe_chol(A, m)
%   [R, QQ_T, info] = isdf.safe_chol(A, m, tol, maxit)
%
%   Only A = X'*X and row count m are required (no explicit X).
%   On convergence QQ_T ≈ I (orthonormal-column Gram matrix); X ≈ Q*R with
%   Q such that Q'*Q = QQ_T, recoverable as Q = X/R when X is available.
%
%   Inputs
%     A     — n-by-n symmetric PSD Gram matrix (X'*X)
%     m     — row count of X (m >= n), used in shift formula
%     tol   — stop when ||QQ_T - I||_F <= tol (default sqrt(n)*eps)
%     maxit — max iterations (default 10)
%
%   Outputs
%     R     — n-by-n upper triangular accumulated factor
%     QQ_T  — current Q'*Q (≈ I when converged)
%     info  — .iters, .fro_err, .shift_count, .x_norm2, .tol

  if nargin < 1 || isempty(A)
    error('isdf:safe_chol:Input', 'A must be a nonempty numeric matrix.');
  end
  if nargin < 2 || isempty(m)
    error('isdf:safe_chol:Input', 'm (number of rows of X) is required.');
  end
  if ~isnumeric(A) || ~ismatrix(A)
    error('isdf:safe_chol:Input', 'A must be a 2-D numeric matrix.');
  end

  if size(A, 1) ~= size(A, 2)
    error('isdf:safe_chol:Size', 'A must be square.');
  end
  n = size(A, 1);
  m = double(m(1));
  if m < n
    error('isdf:safe_chol:Size', 'm must be >= n; got %d < %d.', m, n);
  end

  u = eps(class(A));
  if nargin < 3 || isempty(tol)
    tol = sqrt(n) * u;
  else
    tol = double(tol);
  end
  if nargin < 4 || isempty(maxit)
    maxit = 10;
  else
    maxit = max(1, round(double(maxit)));
  end

  A = double(A);
  % ||X||_2^2 = sigma_max(X'*X) = norm(A, 2) when A = X'*X
  x_norm2_sq = norm(A, 2);
  if ~isfinite(x_norm2_sq) || x_norm2_sq <= 0
    x_norm2_sq = norm(A, 'fro')^2;
  end

  QQ_T = (A + A') / 2;
  R = eye(n);
  shift_count = 0;
  fro_err = inf;

  for it = 1:maxit
    G = (QQ_T + QQ_T') / 2;

    [R_tilde, used_shift] = local_chol_upper(G, m, n, u, x_norm2_sq);
    shift_count = shift_count + double(used_shift);

    % Q_new' * Q_new = R_tilde^{-T} * G * R_tilde^{-1}
    QQ_T = R_tilde' \ (G / R_tilde);
    R = R_tilde * R;

    fro_err = norm(QQ_T - eye(n), 'fro');
    if fro_err <= tol
      break
    end
  end

  if nargout >= 3
    info = struct( ...
      'iters', it, ...
      'fro_err', fro_err, ...
      'shift_count', shift_count, ...
      'x_norm2', sqrt(x_norm2_sq), ...
      'tol', tol);
  end
end

function [R_upper, used_shift] = local_chol_upper(A, m, n, u, x_norm2_sq)
  used_shift = false;
  try
    R_upper = chol(A, 'upper');
  catch
    s = 11 * (m * n + n * (n + 1)) * u * x_norm2_sq;
    if ~isfinite(s) || s <= 0
      s = 11 * (m * n + n * (n + 1)) * u * max(norm(A, 'fro')^2, eps);
    end
    R_upper = chol(A + s * eye(n), 'upper');
    used_shift = true;
  end
end

% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/07/28 ZZ

function varargout = isdf_schur_update(mode, varargin)
%ISDF_SCHUR_UPDATE  Persistent Schur / Cholesky state for adaptive ISDF (double precision).
%
% ---------------------------------------------------------------------------
% Purpose
% ---------------------------------------------------------------------------
% Maintain a growing Gram matrix CCH of ISDF sampling orbitals and its Cholesky
% factorization while adaptiveisdf adds new centroids one-by-one (or in batches).
%
%   CCH ≈ <Psi, Phi> products on the current sampling set (via isdf.prod).
%   L_CCH * L_CCH' = CCH(1:Nisdf,1:Nisdf)   (lower Cholesky)
%   invL_CCH ≈ inv(L_CCH)                   (used for cheap rank-1 Schur updates)
%
% Precision (this package: +adaptive):
%   - Psi_on_grid / Phi_on_grid and isdf.prod I/O: double
%   - CCH / L_CCH / invL_CCH and Schur algebra: double
%   Note: isdf.adaptive.adaptiveisdf currently forwards to +adaptive_double;
%   this file is the double-precision Schur kernel under the +adaptive package.
%   (Contrast +adaptive_single: grids/prod in single, Schur algebra still double.)
%
% ---------------------------------------------------------------------------
% Modes
% ---------------------------------------------------------------------------
%   isdf_schur_update('init', Nisdf, Psi, Phi [, max_cond_number, use_cond_guard])
%       Build CCH = prod(Psi,Phi) on the initial set, Cholesky + invL, seed
%       cch_lambda_max. Preallocates buffers to Nisdfmax = 2*Nisdf.
%
%   ok = isdf_schur_update('update', Nadd, Psi_new, Phi_new [, selected_global_idx])
%       Append Nadd new sampling points. Returns false if the candidate is
%       rejected (non-positive Schur complement or cond guard). See UPDATE
%       section below. Hot path is Nadd == 1.
%
%   [Nisdf, Nadd, invL_CCH, L_CCH, Psi, Phi] = isdf_schur_update('get')
%       Read persistent factors / grids (optional trailing outputs).
%
%   isdf_schur_update('clear')
%       Clear all persistent variables.
%
% ---------------------------------------------------------------------------
% UPDATE (Nadd == 1) — key path used by adaptiveisdf
% ---------------------------------------------------------------------------
% Partition old block (size oldNisdf) and the new row/column:
%
%   CCH_new = [ C11 , C21' ; C21 , C22 ]
%
% where C21 = C2C1H (new vs old products), C22 = C2C2H (new vs new).
% Schur complement of the new diagonal block:
%
%   S = C22 - L21 * L21' ,   L21 = C21 * inv(L11)
%
% For Nadd==1, S is a scalar; L22 = sqrt(S). Then extend
%   L_CCH  and  invL_CCH  with the new row (block inverse of lower-triangular L).
%
% Acceptance checks (any failure => ok=false, state unchanged for that point):
%   1) S > 0  (and L22 finite / positive on MEX paths)
%   2) cond guard: S >= cch_lambda_max / max_cond_number  (if use_cond_guard)
%
% Implementation order for Nadd==1, oldNisdf>0 (first success wins):
%   A) isdf_schur_rank1_prod_mex  — fuse prod + rank-1 Schur in one MEX
%   B) MATLAB: build C2C1H / C2C2H (optionally via MrC1H cache), then
%      isdf_schur_rank1_mex, else pure MATLAB L21 / S / L22 / invL update
%   C) After accept: copy Psi/Phi rows, Nisdf++, refresh cch_lambda_max
%
% Nadd > 1: block Cholesky of S (no cond-guard / MEX fast path); less common.
%           Will be implemented in the future.
%
% Local helpers (below):
%   get_max_diag_LCCH      — Maximum elements of |diag(L_CCH)|
%   isdf_schur_cond_guard_ok — apply max_cond_number gate + stderr message
%

  persistent Nisdf Nadd CCH invL_CCH L_CCH Psi_on_grid Phi_on_grid ...
             Nisdfmax has_rank1_mex has_rank1_prod_mex ...
             cch_lambda_max max_cond_number use_cond_guard

  if nargin < 1
    error('isdf_schur_update:mode', 'First argument ''mode'' is required.');
  end

  switch lower(mode)
    case 'update'
      update_ok = true;
      % varargin{1}: Nadd
      % varargin{2}: Psi on new grid
      % varargin{3}: Phi on new grid
      % varargin{4}: optional global row index for MrC1H cache reuse (Nadd == 1)
      % 
      Nadd = round(double(varargin{1}));
      if (Nadd - varargin{1}) > 1e-6
        error('isdf_schur_update:update', 'Nadd is not an integer.');
      end
      Psi_on_new_grid = double(varargin{2}); % size: Nadd, nb1, nspin
      Phi_on_new_grid = double(varargin{3}); % size: Nadd, nb2, nspin
      selected_global_idx = [];
      if nargin >= 5
        selected_global_idx = varargin{4};
      end
      oldNisdf = Nisdf;
      oldIdx = 1:oldNisdf;
      newIdx = oldNisdf+1:oldNisdf+Nadd;

      if Nadd == 1
        % Fast path: scalar Schur complement update for rank-1 append.
        if oldNisdf > 0
          if isempty(has_rank1_prod_mex)
            has_rank1_prod_mex = ~isempty(which('isdf.adaptive.isdf_schur_rank1_prod_mex'));
          end
          if has_rank1_prod_mex
            try
              [L21_row, L22_mex, invL_new_old] = ...
                isdf.adaptive.isdf_schur_rank1_prod_mex(invL_CCH, oldNisdf, ...
                  Psi_on_new_grid, Phi_on_new_grid, Psi_on_grid, Phi_on_grid);
              if isfinite(L22_mex) && real(L22_mex) > 0 && ...
                 isdf_schur_cond_guard_ok(real(L22_mex)^2, cch_lambda_max, max_cond_number, use_cond_guard, selected_global_idx)
                L_CCH(newIdx, oldIdx) = L21_row;
                L_CCH(newIdx, newIdx) = L22_mex;
                invL22 = double(1.0) / L22_mex;
                invL_CCH(newIdx, newIdx) = invL22;
                invL_CCH(newIdx, oldIdx) = invL_new_old;
                Psi_on_grid(newIdx, :) = Psi_on_new_grid;
                Phi_on_grid(newIdx, :) = Phi_on_new_grid;
                Nisdf = oldNisdf + Nadd;
                cch_lambda_max = get_max_diag_LCCH(invL_CCH, Nisdf, cch_lambda_max);
                if nargout >= 1, varargout{1} = update_ok; end
                return;
              end
            catch
              has_rank1_prod_mex = false;
            end
          end

          if ~isempty(selected_global_idx)
            selected_global_idx = double(selected_global_idx(1));
            wf_data = wave_functions.get();
            isdf.adaptive.MrC1H('ensure', double(wf_data.nc), Nisdfmax);
            cols_have = isdf.adaptive.MrC1H('cached_cols', selected_global_idx);
            if cols_have < oldNisdf
              cstart = cols_have + 1;
              C2C1H_new = isdf.prod(Psi_on_new_grid, Psi_on_grid(cstart:oldNisdf, :), ...
                                    Phi_on_new_grid, Phi_on_grid(cstart:oldNisdf, :));
              isdf.adaptive.MrC1H('set_range', selected_global_idx, cstart, C2C1H_new);
              isdf.adaptive.MrC1H('set_cached_cols', selected_global_idx, oldNisdf);
            end
            C2C1H = isdf.adaptive.MrC1H('get', selected_global_idx, 1, oldNisdf);
          else
            C2C1H = isdf.prod(Psi_on_new_grid, Psi_on_grid(oldIdx, :), Phi_on_new_grid, Phi_on_grid(oldIdx, :));
          end
          C2C2H = isdf.prod(Psi_on_new_grid, Psi_on_new_grid, Phi_on_new_grid, Phi_on_new_grid);

          if isempty(has_rank1_mex)
            has_rank1_mex = ~isempty(which('isdf.adaptive.isdf_schur_rank1_mex'));
          end
          if has_rank1_mex
            try
              [L21_row, L22_mex, invL_new_old] = ...
                isdf.adaptive.isdf_schur_rank1_mex(invL_CCH, oldNisdf, C2C1H, C2C2H);
              if isfinite(L22_mex) && real(L22_mex) > 0 && ...
                 isdf_schur_cond_guard_ok(real(L22_mex)^2, cch_lambda_max, max_cond_number, use_cond_guard, selected_global_idx)
                L_CCH(newIdx, oldIdx) = L21_row;
                L_CCH(newIdx, newIdx) = L22_mex;
                invL22 = double(1.0) / L22_mex;
                invL_CCH(newIdx, newIdx) = invL22;
                invL_CCH(newIdx, oldIdx) = invL_new_old;
                Psi_on_grid(newIdx, :) = Psi_on_new_grid;
                Phi_on_grid(newIdx, :) = Phi_on_new_grid;
                Nisdf = oldNisdf + Nadd;
                cch_lambda_max = get_max_diag_LCCH(invL_CCH, Nisdf, cch_lambda_max);
                if nargout >= 1, varargout{1} = update_ok; end
                return;
              end
            catch
              % Fall back to MATLAB path if mex is unavailable or failed at runtime.
              has_rank1_mex = false;
            end
          end
          invL11 = invL_CCH(oldIdx, oldIdx);
          L21_row = C2C1H * invL11';
        else
          C2C2H = isdf.prod(Psi_on_new_grid, Psi_on_new_grid, Phi_on_new_grid, Phi_on_new_grid);
          L21_row = zeros(1, 0, 'double');
        end
        S = C2C2H - L21_row * L21_row';
        S_real = real(S);
        if S_real <= 0
          if isempty(selected_global_idx)
            fprintf(2, 'isdf_schur_update: non-positive Schur complement S=%.6e (Nadd=1); skip update.\n', S_real);
          else
            fprintf(2, 'isdf_schur_update: non-positive Schur complement S=%.6e at index %d (Nadd=1); skip update.\n', ...
              S_real, round(double(selected_global_idx(1))));
          end
          update_ok = false;
          if nargout >= 1, varargout{1} = update_ok; end
          return;
        end
        if ~isdf_schur_cond_guard_ok(S_real, cch_lambda_max, max_cond_number, use_cond_guard, selected_global_idx)
          update_ok = false;
          if nargout >= 1, varargout{1} = update_ok; end
          return;
        end
        L22 = sqrt(S_real);
        if oldNisdf > 0
          L_CCH(newIdx, oldIdx) = L21_row;
        end
        L_CCH(newIdx, newIdx) = L22;
        invL22 = double(1.0) / L22;
        invL_CCH(newIdx, newIdx) = invL22;
        if oldNisdf > 0
          invL_CCH(newIdx, oldIdx) = -(invL22 * L21_row) * invL11;
        end
      else
        C2C2H = isdf.prod(Psi_on_new_grid, Psi_on_new_grid, Phi_on_new_grid, Phi_on_new_grid);
        if oldNisdf == 0
          C2C1H = zeros(Nadd, 0, 'double');
        else
          C2C1H = isdf.prod(Psi_on_new_grid, Psi_on_grid(oldIdx, :), Phi_on_new_grid, Phi_on_grid(oldIdx, :));
        end
        if oldNisdf > 0
          invL11 = invL_CCH(oldIdx, oldIdx);
          L21 = C2C1H * invL11';
        else
          L21 = zeros(Nadd, 0, 'double');
        end
        L_CCH(newIdx, oldIdx) = L21;
        S = C2C2H - L21 * L21';
        L22 = chol(S, "lower");
        L_CCH(newIdx, newIdx) = L22;
        invL22 = inv(L22);
        invL_CCH(newIdx, newIdx) = invL22;
        if oldNisdf > 0
          invL_CCH(newIdx, oldIdx) = -invL22 * L21 * invL11;
        end
      end

      Psi_on_grid(newIdx, :) = Psi_on_new_grid;
      Phi_on_grid(newIdx, :) = Phi_on_new_grid;
      Nisdf = oldNisdf + Nadd;
      cch_lambda_max = get_max_diag_LCCH(invL_CCH, Nisdf, cch_lambda_max);
      if nargout >= 1, varargout{1} = update_ok; end

    case 'init'
      % varargin{1}: Nisdf
      if nargin < 2 || isempty(varargin{1})
        error('isdf_schur_update:init', '''init'' requires scalar or nonempty Nisdf.');
      end
      Nisdf = double(varargin{1}(1));
      Nisdfmax = Nisdf * 2.0; % TODO: make it a parameter
      %
      nb1 = size(varargin{2}, 2);
      nb2 = size(varargin{3}, 2);
      Psi_on_grid = zeros(Nisdfmax, nb1, 'double');
      Phi_on_grid = zeros(Nisdfmax, nb2, 'double');
      Psi_on_grid(1:Nisdf, :) = double(varargin{2}); % size: Nisdf, nb1
      Phi_on_grid(1:Nisdf, :) = double(varargin{3}); % size: Nisdf, nb2
      
      CCH = zeros(Nisdfmax, Nisdfmax, 'double');
      L_CCH = zeros(Nisdfmax, Nisdfmax, 'double');
      invL_CCH = zeros(Nisdfmax, Nisdfmax, 'double');
      CCH(1:Nisdf, 1:Nisdf) = isdf.prod(Psi_on_grid(1:Nisdf, :), Psi_on_grid(1:Nisdf, :), ...
                      Phi_on_grid(1:Nisdf, :), Phi_on_grid(1:Nisdf, :));
      if condest(CCH(1:Nisdf, 1:Nisdf)) < 1e+12
        L_CCH(1:Nisdf, 1:Nisdf) = chol(CCH(1:Nisdf, 1:Nisdf), "lower");
        invL_CCH(1:Nisdf, 1:Nisdf) = inv(L_CCH(1:Nisdf, 1:Nisdf));
      else
        warning('isdf_schur_update:init', 'CCH is ill-conditioned, using pseudoinverse.');
        L_CCH(1:Nisdf, 1:Nisdf) = chol(CCH(1:Nisdf, 1:Nisdf), "lower");
        invL_CCH(1:Nisdf, 1:Nisdf) = inv(L_CCH(1:Nisdf, 1:Nisdf));
      end
      max_cond_number = 1e12;
      if nargin >= 5 && ~isempty(varargin{4})
        max_cond_number = double(varargin{4});
      end
      if ~isfinite(max_cond_number) || max_cond_number <= 0
        max_cond_number = inf;
      end
      use_cond_guard = true;
      if nargin >= 6 && ~isempty(varargin{5})
        use_cond_guard = logical(varargin{5});
      end
      cch_lambda_max = get_max_diag_LCCH(invL_CCH, Nisdf, 0.0);
      Nadd = Nisdf;
    case 'get'
      varargout{1} = Nisdf;
      varargout{2} = Nadd; 
      varargout{3} = invL_CCH;
      if nargout >= 2
        varargout{4} = L_CCH;
      end
      if nargout >= 5
        varargout{5} = Psi_on_grid;
      end
      if nargout >= 6
        varargout{6} = Phi_on_grid;
      end
      if nargout >= 7
        error('isdf_schur_update:get', 'too much output arguments.');
      end
    case 'clear'
      clear Nisdf Nadd CCH invL_CCH L_CCH Psi_on_grid Phi_on_grid Nisdfmax ...
            has_rank1_mex has_rank1_prod_mex cch_lambda_max max_cond_number use_cond_guard

    otherwise
      error('isdf_schur_update:mode', 'Unknown mode ''%s''.', mode);
  end
end

function lambda_max = get_max_diag_LCCH(invL_CCH, Nisdf, prev_lambda_max)
  if Nisdf <= 0
    lambda_max = prev_lambda_max;
    return;
  end
  d = abs(diag(invL_CCH(1:Nisdf, 1:Nisdf)));
  d = d(isfinite(d) & d > 0);
  if isempty(d)
    lambda_max = prev_lambda_max;
    return;
  end
  lambda_max_now = 1.0 / min(double(d));
  if isempty(prev_lambda_max) || ~isfinite(prev_lambda_max)
    lambda_max = lambda_max_now;
  else
    lambda_max = max(prev_lambda_max, lambda_max_now);
  end
end

function ok = isdf_schur_cond_guard_ok(S_real, lambda_max, max_cond_number, use_cond_guard, selected_global_idx)
  if ~use_cond_guard
    ok = true;
    return;
  end
  if ~isfinite(max_cond_number) || max_cond_number <= 0 || max_cond_number == inf
    ok = true;
    return;
  end
  s_min_allowed = double(lambda_max) / double(max_cond_number);
  ok = isfinite(S_real) && double(S_real) >= s_min_allowed;
  if ~ok
    if isempty(selected_global_idx)
      fprintf(2, 'isdf_schur_update: skip update due to cond guard, S=%.6e < %.6e (lambda_max=%.6e, cond_max=%.6e).\n', ...
        double(S_real), s_min_allowed, double(lambda_max), double(max_cond_number));
    else
      fprintf(2, 'isdf_schur_update: skip index %d due to cond guard, S=%.6e < %.6e (lambda_max=%.6e, cond_max=%.6e).\n', ...
        round(double(selected_global_idx(1))), double(S_real), s_min_allowed, double(lambda_max), double(max_cond_number));
    end
  end
end

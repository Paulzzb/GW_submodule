function varargout = isdf_schur_update(mode, varargin)
%ISDF_SCHUR_UPDATE  Persistent state for Schur / Gram refresh.
%   CCH and invCCH are both Nisdf-by-Nisdf (centroid Gram and its inverse).
%
%   isdf_schur_update('init', Nisdf)
%   isdf_schur_update('init', Nisdf, fftgrid_fine)
%   isdf_schur_update('init', Nisdf, fftgrid_fine, Nremain)
%   [CCH, invCCH, Nremain, fftgrid_fine] = isdf_schur_update('get')
%   isdf_schur_update('clear')
%   ok = isdf_schur_update('update', ...)   % false if Schur complement <= 0 (Nadd==1)

  persistent Nisdf Nadd CCH invL_CCH L_CCH Psi_on_grid Phi_on_grid ...
             Nisdfmax has_rank1_mex has_rank1_prod_mex ...
             cch_lambda_max_proxy max_cond_number use_cond_guard

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
            has_rank1_prod_mex = ~isempty(which('isdf.adaptive_double.isdf_schur_rank1_prod_mex'));
          end
          if has_rank1_prod_mex
            try
              [L21_row, L22_mex, invL_new_old] = ...
                isdf.adaptive_double.isdf_schur_rank1_prod_mex(invL_CCH, oldNisdf, ...
                  Psi_on_new_grid, Phi_on_new_grid, Psi_on_grid, Phi_on_grid);
              if isdf_schur_rank1_ok(L22_mex) && ...
                 isdf_schur_cond_guard_ok(real(L22_mex)^2, cch_lambda_max_proxy, max_cond_number, use_cond_guard, selected_global_idx)
                L_CCH(newIdx, oldIdx) = L21_row;
                L_CCH(newIdx, newIdx) = L22_mex;
                invL22 = double(1.0) / L22_mex;
                invL_CCH(newIdx, newIdx) = invL22;
                invL_CCH(newIdx, oldIdx) = invL_new_old;
                Psi_on_grid(newIdx, :) = Psi_on_new_grid;
                Phi_on_grid(newIdx, :) = Phi_on_new_grid;
                Nisdf = oldNisdf + Nadd;
                cch_lambda_max_proxy = isdf_schur_lambda_max_proxy(invL_CCH, Nisdf, cch_lambda_max_proxy);
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
            isdf.adaptive_double.MrC1H('ensure', double(wf_data.nc), Nisdfmax);
            cols_have = isdf.adaptive_double.MrC1H('cached_cols', selected_global_idx);
            if cols_have < oldNisdf
              cstart = cols_have + 1;
              C2C1H_new = isdf.prod(Psi_on_new_grid, Psi_on_grid(cstart:oldNisdf, :), ...
                                    Phi_on_new_grid, Phi_on_grid(cstart:oldNisdf, :));
              isdf.adaptive_double.MrC1H('set_range', selected_global_idx, cstart, C2C1H_new);
              isdf.adaptive_double.MrC1H('set_cached_cols', selected_global_idx, oldNisdf);
            end
            C2C1H = isdf.adaptive_double.MrC1H('get', selected_global_idx, 1, oldNisdf);
          else
            C2C1H = isdf.prod(Psi_on_new_grid, Psi_on_grid(oldIdx, :), Phi_on_new_grid, Phi_on_grid(oldIdx, :));
          end
          C2C2H = isdf.prod(Psi_on_new_grid, Psi_on_new_grid, Phi_on_new_grid, Phi_on_new_grid);

          if isempty(has_rank1_mex)
            has_rank1_mex = ~isempty(which('isdf.adaptive_double.isdf_schur_rank1_mex'));
          end
          if has_rank1_mex
            try
              [L21_row, L22_mex, invL_new_old] = ...
                isdf.adaptive_double.isdf_schur_rank1_mex(invL_CCH, oldNisdf, C2C1H, C2C2H);
              if isdf_schur_rank1_ok(L22_mex) && ...
                 isdf_schur_cond_guard_ok(real(L22_mex)^2, cch_lambda_max_proxy, max_cond_number, use_cond_guard, selected_global_idx)
                L_CCH(newIdx, oldIdx) = L21_row;
                L_CCH(newIdx, newIdx) = L22_mex;
                invL22 = double(1.0) / L22_mex;
                invL_CCH(newIdx, newIdx) = invL22;
                invL_CCH(newIdx, oldIdx) = invL_new_old;
                Psi_on_grid(newIdx, :) = Psi_on_new_grid;
                Phi_on_grid(newIdx, :) = Phi_on_new_grid;
                Nisdf = oldNisdf + Nadd;
                cch_lambda_max_proxy = isdf_schur_lambda_max_proxy(invL_CCH, Nisdf, cch_lambda_max_proxy);
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
          isdf_schur_update_log_nonpositive(S_real, selected_global_idx);
          update_ok = false;
          if nargout >= 1, varargout{1} = update_ok; end
          return;
        end
        if ~isdf_schur_cond_guard_ok(S_real, cch_lambda_max_proxy, max_cond_number, use_cond_guard, selected_global_idx)
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
      cch_lambda_max_proxy = isdf_schur_lambda_max_proxy(invL_CCH, Nisdf, cch_lambda_max_proxy);
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
      CCH(1:Nisdf, 1:Nisdf) = 0.5*CCH(1:Nisdf, 1:Nisdf) + 0.5*CCH(1:Nisdf, 1:Nisdf)';
      blk = CCH(1:Nisdf, 1:Nisdf);
      [L_try, p_chol] = chol(blk, 'lower');
      if p_chol ~= 0
        jitter = 1e-12 * trace(blk) / max(Nisdf, 1);
        if jitter <= 0
          jitter = 1e-12;
        end
        [L_try, p_chol] = chol(blk + jitter * eye(Nisdf), 'lower');
      end
      if p_chol ~= 0
        error('isdf_schur_update:init:NotPD', ...
          'Initial CCH block (Nisdf=%d) is not positive definite.', Nisdf);
      end
      if condest(blk) < 1e+12
        L_CCH(1:Nisdf, 1:Nisdf) = L_try;
        invL_CCH(1:Nisdf, 1:Nisdf) = inv(L_CCH(1:Nisdf, 1:Nisdf));
      else
        warning('isdf_schur_update:init', 'CCH is ill-conditioned, using pseudoinverse.');
        L_CCH(1:Nisdf, 1:Nisdf) = L_try;
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
      cch_lambda_max_proxy = isdf_schur_lambda_max_proxy(invL_CCH, Nisdf, 0.0);
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
            has_rank1_mex has_rank1_prod_mex cch_lambda_max_proxy max_cond_number use_cond_guard

    otherwise
      error('isdf_schur_update:mode', 'Unknown mode ''%s''.', mode);
  end
end

function ok = isdf_schur_rank1_ok(L22)
  ok = isfinite(L22) && real(L22) > 0;
end

function isdf_schur_update_log_nonpositive(S_real, selected_global_idx)
  if isempty(selected_global_idx)
    fprintf(2, 'isdf_schur_update: non-positive Schur complement S=%.6e (Nadd=1); skip update.\n', S_real);
  else
    fprintf(2, 'isdf_schur_update: non-positive Schur complement S=%.6e at index %d (Nadd=1); skip update.\n', ...
      S_real, round(double(selected_global_idx(1))));
  end
end

function lambda_max_proxy = isdf_schur_lambda_max_proxy(invL_CCH, Nisdf, prev_lambda_max_proxy)
  if Nisdf <= 0
    lambda_max_proxy = prev_lambda_max_proxy;
    return;
  end
  d = abs(diag(invL_CCH(1:Nisdf, 1:Nisdf)));
  d = d(isfinite(d) & d > 0);
  if isempty(d)
    lambda_max_proxy = prev_lambda_max_proxy;
    return;
  end
  lambda_max_now = 1.0 / min(double(d));
  if isempty(prev_lambda_max_proxy) || ~isfinite(prev_lambda_max_proxy)
    lambda_max_proxy = lambda_max_now;
  else
    lambda_max_proxy = max(prev_lambda_max_proxy, lambda_max_now);
  end
end

function ok = isdf_schur_cond_guard_ok(S_real, lambda_max_proxy, max_cond_number, use_cond_guard, selected_global_idx)
  if ~use_cond_guard
    ok = true;
    return;
  end
  if ~isfinite(max_cond_number) || max_cond_number <= 0 || max_cond_number == inf
    ok = true;
    return;
  end
  s_min_allowed = double(lambda_max_proxy) / double(max_cond_number);
  ok = isfinite(S_real) && double(S_real) >= s_min_allowed;
  if ~ok
    if isempty(selected_global_idx)
      fprintf(2, 'isdf_schur_update: skip update due to cond guard, S=%.6e < %.6e (lambda_max_proxy=%.6e, cond_max=%.6e).\n', ...
        double(S_real), s_min_allowed, double(lambda_max_proxy), double(max_cond_number));
    else
      fprintf(2, 'isdf_schur_update: skip index %d due to cond guard, S=%.6e < %.6e (lambda_max_proxy=%.6e, cond_max=%.6e).\n', ...
        round(double(selected_global_idx(1))), double(S_real), s_min_allowed, double(lambda_max_proxy), double(max_cond_number));
    end
  end
end

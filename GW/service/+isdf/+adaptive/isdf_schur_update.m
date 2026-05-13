function varargout = isdf_schur_update(mode, varargin)
%ISDF_SCHUR_UPDATE  Persistent state for Schur / Gram refresh.
%   CCH and invCCH are both Nisdf-by-Nisdf (centroid Gram and its inverse).
%
%   isdf_schur_update('init', Nisdf)
%   isdf_schur_update('init', Nisdf, fftgrid_fine)
%   isdf_schur_update('init', Nisdf, fftgrid_fine, Nremain)
%   [CCH, invCCH, Nremain, fftgrid_fine] = isdf_schur_update('get')
%   isdf_schur_update('clear')

  persistent Nisdf Nadd CCH invL_CCH L_CCH Psi_on_grid Phi_on_grid ...
             Nisdfmax has_rank1_mex has_rank1_prod_mex

  if nargin < 1
    error('isdf_schur_update:mode', 'First argument ''mode'' is required.');
  end

  switch lower(mode)
    case 'update'
      % varargin{1}: Nadd
      % varargin{2}: Psi on new grid
      % varargin{3}: Phi on new grid
      % varargin{4}: optional global row index for MrC1H cache reuse (Nadd == 1)
      % 
      Nadd = round(double(varargin{1}));
      if (Nadd - varargin{1}) > 1e-6
        error('isdf_schur_update:update', 'Nadd is not an integer.');
      end
      Psi_on_new_grid = single(varargin{2}); % size: Nadd, nb1, nspin
      Phi_on_new_grid = single(varargin{3}); % size: Nadd, nb2, nspin
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
              L_CCH(newIdx, oldIdx) = L21_row;
              L_CCH(newIdx, newIdx) = L22_mex;
              invL22 = single(1.0) / L22_mex;
              invL_CCH(newIdx, newIdx) = invL22;
              invL_CCH(newIdx, oldIdx) = invL_new_old;
              Psi_on_grid(newIdx, :) = Psi_on_new_grid;
              Phi_on_grid(newIdx, :) = Phi_on_new_grid;
              Nisdf = oldNisdf + Nadd;
              return;
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
              L_CCH(newIdx, oldIdx) = L21_row;
              L_CCH(newIdx, newIdx) = L22_mex;
              invL22 = single(1.0) / L22_mex;
              invL_CCH(newIdx, newIdx) = invL22;
              invL_CCH(newIdx, oldIdx) = invL_new_old;
              Psi_on_grid(newIdx, :) = Psi_on_new_grid;
              Phi_on_grid(newIdx, :) = Phi_on_new_grid;
              Nisdf = oldNisdf + Nadd;
              return;
            catch
              % Fall back to MATLAB path if mex is unavailable or failed at runtime.
              has_rank1_mex = false;
            end
          end
          invL11 = invL_CCH(oldIdx, oldIdx);
          L21_row = C2C1H * invL11';
          L_CCH(newIdx, oldIdx) = L21_row;
        else
          C2C2H = isdf.prod(Psi_on_new_grid, Psi_on_new_grid, Phi_on_new_grid, Phi_on_new_grid);
          L21_row = zeros(1, 0, 'single');
        end
        S = C2C2H - L21_row * L21_row';
        S_real = real(S);
        if S_real <= 0
          error('isdf_schur_update:update', 'Non-positive Schur complement in Nadd=1 fast path.');
        end
        L22 = sqrt(S_real);
        L_CCH(newIdx, newIdx) = L22;
        invL22 = single(1.0) / L22;
        invL_CCH(newIdx, newIdx) = invL22;
        if oldNisdf > 0
          invL_CCH(newIdx, oldIdx) = -(invL22 * L21_row) * invL11;
        end
      else
        C2C2H = isdf.prod(Psi_on_new_grid, Psi_on_new_grid, Phi_on_new_grid, Phi_on_new_grid);
        if oldNisdf == 0
          C2C1H = zeros(Nadd, 0, 'single');
        else
          C2C1H = isdf.prod(Psi_on_new_grid, Psi_on_grid(oldIdx, :), Phi_on_new_grid, Phi_on_grid(oldIdx, :));
        end
        if oldNisdf > 0
          invL11 = invL_CCH(oldIdx, oldIdx);
          L21 = C2C1H * invL11';
        else
          L21 = zeros(Nadd, 0, 'single');
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
      Psi_on_grid = zeros(Nisdfmax, nb1, 'single');
      Phi_on_grid = zeros(Nisdfmax, nb2, 'single');
      Psi_on_grid(1:Nisdf, :) = single(varargin{2}); % size: Nisdf, nb1
      Phi_on_grid(1:Nisdf, :) = single(varargin{3}); % size: Nisdf, nb2
      
      CCH = zeros(Nisdfmax, Nisdfmax, 'single');
      L_CCH = zeros(Nisdfmax, Nisdfmax, 'single');
      invL_CCH = zeros(Nisdfmax, Nisdfmax, 'single');
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
            has_rank1_mex has_rank1_prod_mex

    otherwise
      error('isdf_schur_update:mode', 'Unknown mode ''%s''.', mode);
  end
end

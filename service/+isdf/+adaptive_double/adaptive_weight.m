function varargout = adaptive_weight(mode, varargin)
%ADAPTIVE_WEIGHT  Persistent state for adaptive weights (skeleton).
%   Mirror of isdf.adaptive_double.isdf_schur_update: init / update / get / clear. Implementation TBD.
%
%   adaptive_weight('init', id)
%   adaptive_weight('update', Nlist [, Nadd])
%   w = adaptive_weight('get')
%   adaptive_weight('set_batch_size', n)
%   adaptive_weight('clear')

  persistent Nw w nbrange1 nbrange2 Phi Psi batch_size

  if nargin < 1
    error('adaptive_weight:mode', 'First argument ''mode'' is required.');
  end

  switch lower(mode)
    case 'set_batch_size'
      if nargin < 2 || isempty(varargin{1})
        error('adaptive_weight:set_batch_size', '''set_batch_size'' requires a positive integer.');
      end
      bs = double(varargin{1});
      if ~(isfinite(bs) && bs >= 1)
        error('adaptive_weight:set_batch_size', 'batch_size must be a positive integer.');
      end
      batch_size = max(1, round(bs));

    case 'init_update'
      % varargin{1}: Nlist; optional varargin{2}: Nadd
      [Nisdf, ~, invL_CCH, ~, Psi_on_grid, Phi_on_grid] ...
              = isdf.adaptive_double.isdf_schur_update('get');
      Nisdfmax = size(Psi_on_grid, 1);
      isdf.adaptive_double.MrC1H('ensure', Nw, Nisdfmax);

      if isempty(batch_size)
        batch_size = 256;
      end
      n_batches = ceil(double(Nw) / double(batch_size));
      for ib = 1:n_batches
        i1 = (ib - 1) * batch_size + 1;
        i2 = min(i1 + batch_size - 1, double(Nw));
        batch_idx = i1:i2;
        batch_idx_col = batch_idx(:);

        cols_have = isdf.adaptive_double.MrC1H('cached_cols', batch_idx_col);
        need_mask = (cols_have < Nisdf);
        if any(need_mask)
          starts = unique(cols_have(need_mask));
          for is = 1:length(starts)
            cstart0 = starts(is);
            row_mask = (cols_have == cstart0);
            rows_now = batch_idx(row_mask);
            cstart = cstart0 + 1;
            MrC1H_new = isdf.prod(Psi(rows_now, :), Psi_on_grid(cstart:Nisdf, :), ...
                                  Phi(rows_now, :), Phi_on_grid(cstart:Nisdf, :));
            isdf.adaptive_double.MrC1H('set_range', rows_now, cstart, MrC1H_new);
            isdf.adaptive_double.MrC1H('set_cached_cols', rows_now, Nisdf);
            cols_have(row_mask) = Nisdf;
          end
        end
        MrC1H = isdf.adaptive_double.MrC1H('get', batch_idx, 1, Nisdf);
        invL1MrC1H = invL_CCH(1:Nisdf, 1:Nisdf) * MrC1H';
        w(batch_idx_col) = w(batch_idx_col) - sum(abs(invL1MrC1H).^2, 1).';
      end

    case 'update'
      % varargin{1}: Nlist; optional varargin{2}: Nadd
      [Nisdf, Nadd, invL_CCH, L_CCH, Psi_on_grid, Phi_on_grid] ...
              = isdf.adaptive_double.isdf_schur_update('get');
      Nisdfmax = size(Psi_on_grid, 1);
      isdf.adaptive_double.MrC1H('ensure', Nw, Nisdfmax);

      if nargin < 2
        Nlist = 1:Nw;
      else
        Nlist = varargin{1};
      end
      if nargin >= 3
        Nadd = round(double(varargin{2}));
      end
        

      Nold = Nisdf - Nadd;
      %
      if Nold == 0
        warning('adaptive_weight:update', 'Nold is 0, which is not expected.');
        warning('use init_update instead');
      end
      %
      if isempty(batch_size)
        batch_size = 256;
      end
      n_list = length(Nlist);
      n_batches = ceil(double(n_list) / double(batch_size));
      for ib = 1:n_batches
        i1 = (ib - 1) * batch_size + 1;
        i2 = min(i1 + batch_size - 1, n_list);
        batch_idx = Nlist(i1:i2);
        batch_idx_col = batch_idx(:);

        if Nold == 0
          MrC1H = isdf.prod(Psi(batch_idx, :), Psi_on_grid(1:Nisdf, :), ...
                            Phi(batch_idx, :), Phi_on_grid(1:Nisdf, :));
          isdf.adaptive_double.MrC1H('set_range', batch_idx, 1, MrC1H);
          isdf.adaptive_double.MrC1H('set_cached_cols', batch_idx_col, Nisdf);
          invL1MrC1H = L_CCH(1:Nisdf, 1:Nisdf) \ MrC1H';
          w(batch_idx_col) = w(batch_idx_col) - sum(abs(invL1MrC1H).^2, 1).';
        else
          cols_have = isdf.adaptive_double.MrC1H('cached_cols', batch_idx_col);
          need_mask = (cols_have < Nold);
          if any(need_mask)
            starts = unique(cols_have(need_mask));
            for is = 1:length(starts)
              cstart0 = starts(is);
              row_mask = (cols_have == cstart0);
              rows_now = batch_idx(row_mask);
              cstart = cstart0 + 1;
              MrC1H_old = isdf.prod(Psi(rows_now, :), Psi_on_grid(cstart:Nold, :), ...
                                    Phi(rows_now, :), Phi_on_grid(cstart:Nold, :));
              isdf.adaptive_double.MrC1H('set_range', rows_now, cstart, MrC1H_old);
              isdf.adaptive_double.MrC1H('set_cached_cols', rows_now, Nold);
              cols_have(row_mask) = Nold;
            end
          end

          MrC2H = isdf.prod(Psi(batch_idx, :), Psi_on_grid(Nold+1:Nisdf, :), ...
                            Phi(batch_idx, :), Phi_on_grid(Nold+1:Nisdf, :));
          isdf.adaptive_double.MrC1H('set_range', batch_idx, Nold + 1, MrC2H);
          isdf.adaptive_double.MrC1H('set_cached_cols', batch_idx_col, Nisdf);
          % invL2MrC2H = L_CCH(Nold+1:Nisdf, Nold+1:Nisdf) \ MrC2H;

          MrC1H = isdf.adaptive_double.MrC1H('get', batch_idx, 1, Nold);
          % invL1MrC1H = L_CCH(1:Nold, 1:Nold) \ MrC1H;

          % L21_invL1MrC1H = L_CCH(Nold+1:Nisdf, 1:Nold) * invL1MrC1H;
          % L21_invL1MrC1H = L_CCH(Nold+1:Nisdf, Nold+1:Nisdf) \ L21_invL1MrC1H;

          % wtilde = sum(abs(invL2MrC2H).^2, 1);
          % wtilde = wtilde - 2 * real(sum(conj(L21_invL1MrC1H) .* invL2MrC2H, 1));
          % wtilde = wtilde + sum(abs(L21_invL1MrC1H).^2, 1);

          % invL2MrC2H = double(invL_CCH(Nold+1:Nisdf, Nold+1:Nisdf)) * double(MrC2H);
          % invL21MrC1H = double(invL_CCH(Nold+1:Nisdf, 1:Nold)) * double(MrC1H);
          % wtilde2 = double(sum(abs(invL2MrC2H).^2, 1));
          % wtilde2 = wtilde2 + double(2 * real(sum(conj(invL21MrC1H) .* invL2MrC2H, 1)));
          % wtilde2 = wtilde2 + double(sum(abs(invL21MrC1H).^2, 1));

          invL2MrC2H = invL_CCH(Nold+1:Nisdf, Nold+1:Nisdf) * MrC2H';
          invL21MrC1H = invL_CCH(Nold+1:Nisdf, 1:Nold) * MrC1H';
          wtilde3 = sum(abs(invL2MrC2H + invL21MrC1H).^2, 1);
          % wtilde4 = zeros(size(wtilde3));
          % for i = 1:length(wtilde3)
          %   wtilde4(i) = norm(invL2MrC2H(:, i) + invL21MrC1H(:, i))^2;
          % end
          % if (norm(wtilde3 - wtilde4) / norm(wtilde3) > 1e-3)
          %   error('adaptive_weight:update', 'Weight update is not consistent. Output: %f', norm(wtilde2 - wtilde3) / norm(wtilde2));
          % end
          % invL2MrC2H = double(invL_CCH(Nold+1:Nisdf, Nold+1:Nisdf)) * double(MrC2H);
          % invL21MrC1H = double(invL_CCH(Nold+1:Nisdf, 1:Nold)) * double(MrC1H);
          % wtilde2 = sum(abs(invL2MrC2H).^2, 1);
          % wtilde2 = wtilde2 + 2 * real(sum(conj(invL21MrC1H) .* invL2MrC2H, 1));
          % wtilde2 = wtilde2 + sum(abs(invL21MrC1H).^2, 1);
          % output = norm(wtilde - double(wtilde2)) / norm(wtilde);
          % if output > 1e-5
          %   error('adaptive_weight:update', 'Weight update is not consistent. Output: %f', output);
          % end
          % w(batch_idx_col) = w(batch_idx_col) - wtilde.';
          w(batch_idx_col) = w(batch_idx_col) - wtilde3.';
          if any(wtilde3 < 0)
            warning('adaptive_weight:update', 'Weight is negative. This is not expected.');
          end
        end
      end
      % diff = norm( w(Nlist) );
      % if diff > 1e-5
      %   error('adaptive_weight:update', 'When updating the weight, covered Nlist value should be 0.0.');
      % end

    case 'init'
      if nargin < 2 || isempty(varargin{1})
        error('adaptive_weight:init', '''init'' requires isdf id as second argument.');
      end
      id = varargin{1};
      wf_data = wave_functions.get();
      isdf_data = isdf.get(id);
      if isempty(isdf_data.nrange1) || isempty(isdf_data.nrange2)
        error('adaptive_weight:init:nrange', ...
          ['Missing cached nrange in ISDF id=%d. ' ...
           'Build it first via isdf.set_nrange(id, config.SYSTEM).'], ...
          int32(id));
      end
      Nw = double(wf_data.nc);
      w = zeros(Nw, 1, 'double');
      nbrange1 = double(isdf_data.nrange1);
      nbrange2 = double(isdf_data.nrange2);

      k_data = lattice.manager('k', 'get');
      Psi = zeros(Nw, numel(nbrange1) * double(k_data.nbz), 'double');
      Phi = zeros(Nw, numel(nbrange2) * double(k_data.nbz), 'double');

      count1 = 1;
      count2 = 1;
      ispin = 1;
      for ikbz = 1:k_data.nbz
        ikibz = k_data.bz2ibz(ikbz, 1);
        ikrot = k_data.bz2rot(ikbz, 1);
        for ib1 = 1:numel(nbrange1)
          param1 = [nbrange1(ib1), ikibz, ikrot, ispin];
          Psi(:, count1) = double(wave_functions.WF_apply_symm(param1));
          count1 = count1 + 1;
        end
        for ib2 = 1:numel(nbrange2)
          param2 = [nbrange2(ib2), ikibz, ikrot, ispin];
          Phi(:, count2) = double(wave_functions.WF_apply_symm(param2));
          count2 = count2 + 1;
        end
      end

      isdf.adaptive_double.MrC1H('clear');
      % Calculate diag(MMH), fill them into w
      for ic = 1:Nw
        Psi_r = Psi(ic, :);
        Phi_r = Phi(ic, :);
        w(ic) = real(isdf.prod(Psi_r, Psi_r, Phi_r, Phi_r));
      end
      

    case 'get'
      % varargout: TBD (e.g. weight vector)
      Nlist = 1:Nw;
      if nargin >= 2
        Nlist = varargin{1};
      end
      varargout{1} = w(Nlist);

    case 'clear'
      isdf.adaptive_double.MrC1H('clear');
      clear Nw w nbrange1 nbrange2 Phi Psi batch_size

    case 'get_wf_xga'
      if isempty(Psi) || isempty(Phi)
        error('adaptive_weight:get_wf_xga', ...
          'Psi/Phi not initialized. Call adaptive_weight(''init'', id) first.');
      end
      Nlist = varargin{1};
      varargout{1} = Psi(Nlist, :);
      varargout{2} = Phi(Nlist, :);

    otherwise
      error('adaptive_weight:mode', 'Unknown mode ''%s''.', mode);
  end
end

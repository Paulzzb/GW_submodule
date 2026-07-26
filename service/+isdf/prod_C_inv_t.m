function varargout = prod_C_inv_t(action, varargin)
%PROD_C_INV_T  Apply truncated V*Lambda^{-ratio} factors (left/right).
%
%   l_keep = isdf.prod_C_inv_t('set', CCHq)
%   l_keep = isdf.prod_C_inv_t('set', CCHq, trunc)
%   l_keep = isdf.prod_C_inv_t('set', CCHq, trunc, ratio)
%   isdf.prod_C_inv_t('set_factors', V_trunc, Lambda_trunc, ratio)
%   [V_trunc, Lambda_trunc, l_keep, ratio] = isdf.prod_C_inv_t('get_factors')
%   Y = isdf.prod_C_inv_t('prod', 'r', A)   % A * V_trunc * Lambda_trunc^{ratio-1} (padded to N columns)
%   Y = isdf.prod_C_inv_t('prod', 'l', A)   % Lambda_trunc^{-ratio} * V_trunc' * A (padded to N rows)
%   [L_full, R_full] = isdf.prod_C_inv_t('assemble')
%   isdf.prod_C_inv_t('clear')
%
% Truncation rule keeps sigma_i > trunc * sigma_1.

  persistent V_trunc lambda_trunc lam_pow_l lam_pow_r l_keep Ndim ratio_used

  if nargin < 1
    error('prod_C_inv_t:action', 'First argument ''action'' is required.');
  end

  switch lower(action)
    case 'set'
      if nargin < 2 || isempty(varargin{1})
        error('prod_C_inv_t:set', '''set'' requires CCHq.');
      end
      CCHq = varargin{1};
      trunc = 0;
      ratio = 0.5;
      if nargin >= 3 && ~isempty(varargin{2})
        trunc = double(varargin{2});
      end
      if nargin >= 4 && ~isempty(varargin{3})
        ratio = double(varargin{3});
      end
      ratio_used = ratio;

      CCHq = (CCHq + CCHq') / 2;
      [V, sigma] = eig(CCHq, 'vector');
      [sigma, perm] = sort(sigma, 'descend');
      V = V(:, perm);
      Ndim = size(CCHq, 1);

      if isempty(sigma) || sigma(1) <= 0
        l_keep = 0;
        V_trunc = zeros(Ndim, 0);
        lambda_trunc = zeros(0, 1);
        lam_pow_l = zeros(0, 1);
        lam_pow_r = zeros(0, 1);
      else
        keep = find(sigma > sigma(1) * trunc);
        l_keep = numel(keep);
        V_trunc = V(:, keep);
        lambda_trunc = sigma(keep);
        lam_pow_l = lambda_trunc .^ (ratio-1);
        lam_pow_r = lambda_trunc .^ (-ratio);
      end

      if nargout >= 1
        varargout{1} = l_keep;
      end

    case 'set_factors'
      if nargin < 3
        error('prod_C_inv_t:set_factors', ...
          '''set_factors'' requires V_trunc and Lambda_trunc.');
      end
      V_trunc = varargin{1};
      lambda_trunc = varargin{2};
      ratio = 0.5;
      if nargin >= 4 && ~isempty(varargin{3})
        ratio = double(varargin{3});
      end
      ratio_used = ratio;
      Ndim = size(V_trunc, 1);
      l_keep = size(V_trunc, 2);
      if numel(lambda_trunc) ~= l_keep
        error('prod_C_inv_t:set_factors', ...
          'Lambda_trunc length (%d) must match size(V_trunc,2) (%d).', numel(lambda_trunc), l_keep);
      end
      lambda_trunc = lambda_trunc(:);
      lam_pow_l = lambda_trunc .^ (ratio-1);
      lam_pow_r = lambda_trunc .^ (-ratio);

    case 'get_factors'
      if isempty(V_trunc)
        error('prod_C_inv_t:get_factors', ...
          'Factors not set. Call ''set'' or ''set_factors'' first.');
      end
      varargout{1} = V_trunc;
      if nargout >= 2
        varargout{2} = lambda_trunc;
      end
      if nargout >= 3
        varargout{3} = l_keep;
      end
      if nargout >= 4
        varargout{4} = ratio_used;
      end

    case 'prod'
      if isempty(V_trunc)
        error('prod_C_inv_t:prod', ...
          'Factors not set. Call prod_C_inv_t(''set'', CCHq, trunc, ratio) first.');
      end
      if nargin < 3
        error('prod_C_inv_t:prod', '''prod'' requires side (''l'' or ''r'') and matrix A.');
      end
      side = lower(string(varargin{1}));
      A = varargin{2};

      if side == "r"
        % A * V_trunc * Lambda_trunc^{ratio-1}, padded to Ndim columns.
        Y = zeros(size(A, 1), Ndim);
        if l_keep > 0
          Y(:, 1:l_keep) = (A * V_trunc) .* (lam_pow_r.');
        end
        varargout{1} = Y;
      elseif side == "l"
        % Lambda_trunc^{-ratio} * V_trunc' * A, padded to Ndim rows.
        Y = zeros(Ndim, size(A, 2));
        if l_keep > 0
          Y(1:l_keep, :) = (lam_pow_l .* (V_trunc' * A));
        end
        varargout{1} = Y;
      else
        error('prod_C_inv_t:prod', 'side must be ''l'' or ''r'' (got ''%s'').', side);
      end

    case 'assemble'
      if isempty(V_trunc)
        error('prod_C_inv_t:assemble', ...
          'Factors not set. Call prod_C_inv_t(''set'', CCHq, trunc, ratio) first.');
      end
      L = zeros(Ndim, Ndim);
      R = zeros(Ndim, Ndim);
      if l_keep > 0
        L(1:l_keep, :) = (lam_pow_l .* V_trunc');
        R(:, 1:l_keep) = V_trunc .* (lam_pow_r.');
      end
      if nargout >= 1
        varargout{1} = L;
      end
      if nargout >= 2
        varargout{2} = R;
      end

    case 'clear'
      clear V_trunc lambda_trunc lam_pow_l lam_pow_r l_keep Ndim ratio_used

    otherwise
      error('prod_C_inv_t:action', 'Unknown action ''%s''.', action);
  end
end

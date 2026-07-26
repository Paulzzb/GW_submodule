function varargout = prod_C_inv_sqrt(action, varargin)
%PROD_C_INV_SQRT  Backward-compatible wrapper to prod_C_inv_t (ratio=1/2).
%
% Prefer using isdftest.prod_C_inv_t directly.

  switch lower(action)
    case 'set'
      CCHq = varargin{1};
      s_cut = 0;
      if nargin >= 3 && ~isempty(varargin{2})
        s_cut = varargin{2};
      end
      [varargout{1:nargout}] = isdftest.prod_C_inv_t('set', CCHq, s_cut, 0.5);

    case 'prod'
      [varargout{1:nargout}] = isdftest.prod_C_inv_t('prod', varargin{:});

    case 'assemble'
      % Legacy API returned a full N-by-N operator. Reconstruct from L/R factors.
      [L, R] = isdftest.prod_C_inv_t('assemble');
      varargout{1} = R * L;

    case 'clear'
      isdftest.prod_C_inv_t('clear');

    otherwise
      error('prod_C_inv_sqrt:action', 'Unknown action ''%s''.', action);
  end
end

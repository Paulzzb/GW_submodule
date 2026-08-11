%
% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08
function coeff = isdf_get_coeff(id, c1, c2, varargin)
    persistent firsttime is_coarse N_MAX

  if nargin == 1 && (ischar(id) || (isstring(id) && isscalar(id)))
    cmd = lower(string(id));
    if cmd == "reset"
      firsttime = [];
      is_coarse = [];
      coeff = [];
      return;
    end
  end

  if isempty(firsttime)
    N_MAX = isdf.manager('nmax');
    firsttime = false(N_MAX, 1);
    is_coarse = false(N_MAX, 1);
  end

  if ~firsttime(id)
    isdf_data = isdf.get(id);
    type = isdf_data.interp_scheme;
    if strcmp(type, 'coarse'); is_coarse = true; end
    firsttime(id) = false;
  end
  
  if is_coarse(id)
    % 鍒ゆ柇varargin鏄惁涓虹┖锛堝簲鏈変袱涓€硷級
    if numel(varargin) ~= 2
      error('isdf_get_coeff:InvalidInput', 'Expected two arguments for coarse interpolation');
    end
    irot = varargin{1};
    jrot = varargin{2};
    if irot == 1; c1 = c1; else; c1 = isdf.coeff.isdf_apply_symm_on_coarse(id, c1, irot); end
    if jrot == 1; c2 = c2; else; c2 = isdf.coeff.isdf_apply_symm_on_coarse(id, c2, jrot); end
    coeff = conj(c1) .* c2; 
  else
    coeff = conj(c1) .* c2;
  end
end



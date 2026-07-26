% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/06/02

function d_out = convert_isdf_m_to_isdftest_m(d_in)
%CONVERT_ISDF_M_TO_ISDFTEST_M  Copy shared fields from +isdf pool object to +isdftest.
%
%   d_out = isdftest.adaptive.convert_isdf_m_to_isdftest_m(d_in)
%
% isdftest-only fields (CCHq, CCHq_inv_sqrt, svd_s_cut) keep class defaults.

  if ~isa(d_in, 'isdf.base.isdf_m')
    error('isdftest:adaptive:convertType', ...
      'Input must be isdf.base.isdf_m (got %s).', class(d_in));
  end

  d_out = isdftest.base.isdftest_m();
  src_props = properties(d_in);
  for k = 1:numel(src_props)
    p = src_props{k};
    if isprop(d_out, p)
      d_out.(p) = d_in.(p);
    end
  end
end

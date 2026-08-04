% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/04 ZZ

function numerical_cond_report(action, varargin)
%NUMERICAL_COND_REPORT  Compatibility shim -> isdf.report.cond.
%
%   Prefer: isdf.report.cond(...)
%   See +report/NAMING.md.

  isdf.report.cond(action, varargin{:});
end

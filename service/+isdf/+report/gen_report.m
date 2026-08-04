% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/04 ZZ

function fpath = gen_report(report, outDir)
%GEN_REPORT  Compatibility shim -> isdf.report.hf.
%
%   Prefer: isdf.report.hf(...)
%   See +report/NAMING.md.

  if nargin < 2
    fpath = isdf.report.hf(report);
  else
    fpath = isdf.report.hf(report, outDir);
  end
end

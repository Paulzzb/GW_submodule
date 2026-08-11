% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/07 ZZ

function p = get_report_path()
%GET_REPORT_PATH  Absolute path of the open report file, or ''.
  s = state_('get');
  p = '';
  if isstruct(s) && isfield(s, 'report_path')
    p = char(string(s.report_path));
  end
end

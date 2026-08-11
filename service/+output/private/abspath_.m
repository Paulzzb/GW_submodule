% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/07 ZZ

function p = abspath_(path)
%ABSPATH_  Resolve a possibly-relative path against pwd.
%
%   Relative report/log paths are otherwise ambiguous across cd() and
%   across ensure_init no-ops that compare path strings only.

  path = char(string(path));
  if isempty(path)
    p = '';
    return
  end

  if path(1) == filesep
    p = path;
    return
  end
  if ispc && numel(path) >= 2 && path(2) == ':'
    p = path;
    return
  end

  p = fullfile(pwd, path);
end

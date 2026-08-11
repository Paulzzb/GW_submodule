% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/01 ZZ

function close(name)
%CLOSE  Close a named output file opened by output.open.
%
%   output.close('qp')
%   output.close()          % close all named files

  s = state_('get');

  if nargin < 1 || isempty(name)
    for i = 1:numel(s.named)
      if s.named(i).fid > 0
        try, fclose(s.named(i).fid); catch, end %#ok<CTCH>
      end
    end
    s.named = struct('name', {}, 'path', {}, 'fid', {});
    state_('set', s);
    return
  end

  name = char(string(name));
  keep = true(1, numel(s.named));
  for i = 1:numel(s.named)
    if strcmp(s.named(i).name, name)
      if s.named(i).fid > 0
        try, fclose(s.named(i).fid); catch, end %#ok<CTCH>
      end
      keep(i) = false;
    end
  end
  s.named = s.named(keep);
  state_('set', s);
end

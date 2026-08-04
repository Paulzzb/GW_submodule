% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/01 ZZ

function open(name, path, mode)
%OPEN  Register and open a named output file (Yambo OF_open_close).
%
%   output.open('qp', '<case>/o.qp')
%   output.open('qp', '<case>/o.qp', 'a')   % append
%
%   If path already exists and mode is not append, the existing file is
%   rotated to path_01, path_02, ... before opening a fresh file.
%   Subsequent output.msg('o qp', ...) writes to this file.

  if nargin < 2
    error('output:open:Nargin', 'output.open(name, path) requires both args.');
  end
  if nargin < 3 || isempty(mode)
    mode = 'w';
  end

  name = char(string(name));
  path = char(string(path));
  mode = char(string(mode));

  if isempty(name) || isempty(path)
    error('output:open:Empty', 'name and path must be non-empty.');
  end

  s = state_('get');

  % Already open under this name: close first
  for i = 1:numel(s.named)
    if strcmp(s.named(i).name, name)
      if s.named(i).fid > 0
        try, fclose(s.named(i).fid); catch, end %#ok<CTCH>
      end
      s.named(i) = [];
      break
    end
  end

  d = fileparts(path);
  if ~isempty(d) && ~isfolder(d)
    mkdir(d);
  end

  is_append = any(mode == 'a');
  if ~is_append && isfile(path)
    local_rename_if_exists(path);
  end

  if is_append
    fopen_mode = 'a';
  else
    fopen_mode = 'w';
  end

  [fid, errmsg] = fopen(path, fopen_mode);
  if fid < 0
    error('output:open:OpenFailed', 'Cannot open %s: %s', path, errmsg);
  end

  entry.name = name;
  entry.path = path;
  entry.fid = fid;
  s.named(end+1) = entry; %#ok<AGROW>
  state_('set', s);
end

function local_rename_if_exists(path)
  [folder, base, ext] = fileparts(path);
  for k = 1:99
    cand = fullfile(folder, sprintf('%s_%02d%s', base, k, ext));
    if ~isfile(cand)
      movefile(path, cand);
      return
    end
  end
end

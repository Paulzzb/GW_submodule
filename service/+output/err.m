% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/01 ZZ

function err(text, varargin)
%ERR  Emit a formatted error (screen + report + log) and halt.
%
%   output.err('missing groundstate')
%   output.err('bad value %d', n)

  if nargin < 1
    text = '';
  end

  if nargin < 2 || isempty(varargin)
    body = char(string(text));
  else
    body = sprintf(char(string(text)), varargin{:});
  end

  st = dbstack;
  loc = 'in anonymous context';
  skip = {'err', 'output.err'};
  for i = 1:numel(st)
    if ~any(strcmp(st(i).name, skip))
      loc = sprintf('in %s at line %d', st(i).name, st(i).line);
      break
    end
  end

  line1 = sprintf('[QP-ERROR] %s', body);
  line2 = sprintf('         --> %s', loc);

  % Report/log via emit; stderr always for interactive visibility.
  flags = parse_how_('v0r');
  flags.log = true;
  emit_(flags, line1);
  emit_(flags, line2);
  fprintf(2, '%s\n%s\n', line1, line2);

  error('[QP] Execution terminated.');
end

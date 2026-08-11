% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/01 ZZ

function line = format_line_(text, args)
%FORMAT_LINE_  Build a display line with optional sprintf args and tag/indent.

  s = state_('get');

  if nargin < 2 || isempty(args)
    body = char(string(text));
  else
    body = sprintf(char(string(text)), args{:});
  end

  prefix = '';
  if s.showtag && ~isempty(s.module_stack)
    prefix = sprintf('[%s] ', s.module_stack{end});
  end

  indent = '';
  if s.depth >= 0
    indent = repmat(' ', 1, s.depth);
  end

  line = [indent, prefix, body];
end

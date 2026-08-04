% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/01 ZZ

function warn(text, varargin)
%WARN  Emit a warning to screen + report + log (always, verbose-gated at 0).
%
%   output.warn('units of wfncut not checked')
%   output.warn('missing field %s', name)

  if nargin < 1
    text = '';
  end
  line = format_line_(['[WARN] ', char(string(text))], varargin);
  flags = parse_how_('v0rs');
  flags.log = true;
  emit_(flags, line);
end

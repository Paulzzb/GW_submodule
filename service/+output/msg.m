% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/01 ZZ

function msg(how, text, varargin)
%MSG  Unified message writer (Yambo-like how destinations).
%
%   output.msg(how, text)
%   output.msg(how, fmt, A, B, ...)
%
%   how: combination of
%     s  screen
%     r  report file
%     l  log file
%     o <name>  named output file (opened via output.open)
%     n  blank line before (prefix) / after (suffix)
%     vN verbosity gate N=0|1|2
%
%   Examples:
%     output.msg('rs', 'Exchange Self-Energy')
%     output.msg('r',  'ISDF ratio : %.2f', x)
%     output.msg('v2l', 'iq=%d iter=%d', iq, it)
%     output.msg('o qp', '%4d %12.6f', n, e)

  if nargin < 1
    error('output:msg:Nargin', 'output.msg(how, text, ...) requires how.');
  end
  if nargin < 2
    text = '';
  end

  flags = parse_how_(how);
  line = format_line_(text, varargin);
  emit_(flags, line);
end

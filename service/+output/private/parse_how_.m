% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/01 ZZ

function flags = parse_how_(how)
%PARSE_HOW_  Parse a Yambo-like how string into destination flags.
%
%   how characters:
%     s          screen
%     r          report file
%     l          log file
%     o <name>   named output file (rest of string after leading 'o')
%     n          blank line before (leading n) / after (trailing n)
%     v<N>       verbosity gate N = 0,1,2 (default 0 = always)
%
%   Combinations: 'rs', 'nr', 'rn', 'nrs', 'v2l', 'o qp', ...

  flags = struct( ...
    'screen', false, ...
    'report', false, ...
    'log', false, ...
    'of_name', '', ...
    'blank_before', false, ...
    'blank_after', false, ...
    'verbose', 0 ...
    );

  if nargin < 1 || isempty(how)
    return
  end

  how = char(string(how));
  how = strtrim(how);
  if isempty(how)
    return
  end

  % Named output: "o <name>" or "o<name>"
  if how(1) == 'o'
    rest = strtrim(how(2:end));
    if isempty(rest)
      error('output:parse_how:MissingOfName', ...
        'how ''o'' requires a file name, e.g. ''o qp''.');
    end
    flags.of_name = rest;
    return
  end

  % Verbosity gate: v0 / v1 / v2 anywhere in the control token
  tok = how;
  sp = find(how == ' ', 1, 'first');
  if ~isempty(sp)
    tok = how(1:sp-1);
  end

  vm = regexp(tok, 'v([0-2])', 'tokens', 'once');
  if ~isempty(vm)
    flags.verbose = str2double(vm{1});
    tok = regexprep(tok, 'v[0-2]', '');
  end

  if ~isempty(tok) && tok(1) == 'n'
    flags.blank_before = true;
    tok = tok(2:end);
  end
  if ~isempty(tok) && tok(end) == 'n'
    flags.blank_after = true;
    tok = tok(1:end-1);
  end

  flags.report = any(tok == 'r');
  flags.log = any(tok == 'l');
  flags.screen = any(tok == 's');
end

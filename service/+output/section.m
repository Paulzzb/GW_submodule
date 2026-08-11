% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/01 ZZ

function section(mode, name)
%SECTION  Print a numbered section header (Yambo COM_section style).
%
%   output.section('*', 'Groundstate loading')   % top-level [01]
%   output.section('+', 'ISDF coefficients')     % sub-level [01.02]
%   output.section('-')                          % pop one level (+elapsed)
%   output.section('r')                          % reset depth
%
%   Screen policy:
%     '*' headers always go to screen+report (how 'nrs')
%     '+'/'-' headers go to report; screen only if verbose >= 2

  if nargin < 1 || isempty(mode)
    error('output:section:Nargin', 'output.section(mode, name) requires mode.');
  end
  if nargin < 2
    name = '';
  end

  mode = char(string(mode));
  name = char(string(name));
  s = state_('get');

  switch mode
    case 'r'
      s.depth = -1;
      s.isec = zeros(1, 5);
      s.sec_tics = {};
      state_('set', s);
      return

    case '*'
      if s.depth >= 0 && ~isempty(s.sec_tics)
        local_emit_elapsed(s, 0);
        s = state_('get');
      end
      s.depth = 0;
      s.isec(2:end) = 0;
      s.isec(1) = s.isec(1) + 1;
      s.sec_tics = {tic};

    case '+'
      s.depth = s.depth + 1;
      if s.depth > 4
        s.depth = 4;
      end
      s.isec(s.depth+1) = s.isec(s.depth+1) + 1;
      if s.depth+1 < numel(s.isec)
        s.isec(s.depth+2:end) = 0;
      end
      s.sec_tics{s.depth+1} = tic; %#ok<AGROW>

    case '-'
      if s.depth >= 0
        local_emit_elapsed(s, s.depth);
        s = state_('get');
        s.isec(s.depth+1:end) = 0;
        if numel(s.sec_tics) >= s.depth+1
          s.sec_tics = s.sec_tics(1:s.depth);
        end
        s.depth = s.depth - 1;
      end
      state_('set', s);
      return

    otherwise
      error('output:section:BadMode', ...
        'Unknown section mode ''%s'' (use *, +, -, r).', mode);
  end

  state_('set', s);

  if isempty(strtrim(name))
    return
  end

  secnm = local_secnm(s.isec);
  header = sprintf('%s %s', secnm, name);
  rule = repmat('=', 1, max(8, min(72, numel(header))));

  how_hdr = local_how_for_mode(mode, s.verbose);
  output.msg(how_hdr, '%s', header);
  output.msg(regexprep(how_hdr, 'n', ''), '%s', rule);
end

function secnm = local_secnm(isec)
  secnm = sprintf('[%02d', isec(1));
  for i = 2:numel(isec)
    if isec(i) == 0
      break
    end
    secnm = sprintf('%s.%02d', secnm, isec(i));
  end
  secnm = [secnm, ']'];
end

function how = local_how_for_mode(mode, verbose)
  switch mode
    case '*'
      how = 'nrs';
    case '+'
      if verbose >= 2
        how = 'nrs';
      else
        how = 'nr';
      end
    otherwise
      how = 'r';
  end
end

function local_emit_elapsed(s, depth)
  if depth < 0 || numel(s.sec_tics) < depth+1 || isempty(s.sec_tics{depth+1})
    return
  end
  elapsed = toc(s.sec_tics{depth+1});
  how = 'r';
  if s.verbose >= 2
    how = 'rs';
  end
  flags = parse_how_(how);
  line = format_line_(sprintf('Timing: %.3f s', elapsed), {});
  emit_(flags, line);
end

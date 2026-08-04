% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08
%
% Port of Yambo TIMING_string.F: format CPU seconds as d/h/m/s pieces.
% MATLAB R2008a: char output only.

function s = timing_string(tcpu)
  ltcpu = abs(double(tcpu));
  d = floor(ltcpu / 86400);
  ltcpu = ltcpu - d * 86400;
  h = floor(ltcpu / 3600);
  ltcpu = ltcpu - h * 3600;
  m = floor(ltcpu / 60);
  ssec = floor(ltcpu - m * 60);

  cd = '';
  ch = '';
  cm = '';
  cs = '';

  if d > 0
    cd = sprintf('%02dd', d);
    ch = sprintf('-%02dh', h);
    cm = sprintf('-%02dm', m);
  elseif h > 0
    ch = sprintf('%02dh', h);
    cm = sprintf('-%02dm', m);
  elseif m > 0
    cm = sprintf('%02dm', m);
    cs = sprintf('-%02ds', ssec);
  elseif ssec > 1
    cs = sprintf('%02ds', ssec);
  else
    cs = sprintf('%ds', max(0, ssec));
  end

  s = [cd, ch, cm, cs];
end

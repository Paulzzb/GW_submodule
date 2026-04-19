% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08
%
% Port of Yambo LIVE_timing_update.F (TTY path: one line, hash bar).

function tm = live_timing_update(tm, xch, ech)
  liv = tm.live;

  hn = double(liv.hashes_now);
  hd = double(liv.hashes_done);
  sd = double(liv.steps_done);
  ts = double(liv.time_steps);
  nh = double(liv.nhash);

  if hn == hd && sd ~= 0
    tm.live = liv;
    return;
  end

  nm = liv.timing_name;
  if isempty(strtrim(char(nm)))
    tm.live = liv;
    return;
  end

  if hn ~= hd
    liv.hashes_done = int32(hn);
    hd = double(liv.hashes_done);
  end

  if ts > 0
    perc = floor(sd / ts * 100);
  else
    perc = 0;
  end
  if perc > 100
    perc = 100;
  end

  bar_hashes = repmat('#', 1, hd);
  nsp = max(0, nh - hd);
  bar_spaces = repmat(' ', 1, nsp);

  line = sprintf('%s |%s%s| [%03d%%] %s(E) %s(X)', char(nm), bar_hashes, bar_spaces, perc, char(ech), char(xch));
  fprintf('%s\n', line);

  tm.live = liv;
end

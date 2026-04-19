% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08
%
% Port of Yambo LIVE_timing_close.F (flush omitted: no log unit).

function tm = live_timing_close(tm)
  liv = tm.live;
  if ~liv.live_timing_is_on
    return;
  end

  ts = double(liv.time_steps);
  sd = double(liv.steps_done);
  if ts > 0 && sd ~= ts
    tm = live_timing_add(tm, ts - sd);
    liv = tm.live;
  end

  liv.live_timing_is_on = false;
  tm.live = liv;
end

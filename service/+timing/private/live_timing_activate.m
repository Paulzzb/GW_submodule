% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08
%
% Port of Yambo LIVE_timing_activate.F

function tm = live_timing_activate(tm, name, steps, depth)
  liv = tm.live;

  liv.live_timing_is_on = true;
  liv.timing_name = strtrim(char(name));
  liv.hashes_done = int32(0);
  liv.hashes_now = int32(0);
  liv.time_steps = int32(steps);
  liv.steps_done = int32(0);
  liv.steps_done_in_the_memory = int32(0);
  liv.cput_last_report = 0;
  liv.cput_last_estimate = 0;
  liv.memory_steps = int32(0);

  if nargin >= 4 && ~isempty(depth)
    if isnumeric(depth) && isfinite(depth)
      liv.memory_steps = int32(floor(double(steps) * double(depth)));
    end
  end

  tm.live = liv;
  tm = timing.timing_get_time(tm, 'INIT_SEG');
  tm = live_timing_update(tm, '--', '--');
end

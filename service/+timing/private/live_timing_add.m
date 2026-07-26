% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08
%
% Port of Yambo LIVE_timing_add.F (serial MATLAB).

function tm = live_timing_add(tm, steps)
  liv = tm.live;
  if ~liv.live_timing_is_on
    return;
  end

  steps = double(steps);
  tm = timing.timing_get_time(tm, 'SEG');
  liv = tm.live;

  rts = double(liv.live_report_min_seconds);
  if rts < 0
    rts = 0;
  end
  time_steps = double(liv.time_steps);
  steps_done = double(liv.steps_done);

  if time_steps <= 0
    tm.live = liv;
    return;
  end

  if steps_done + steps <= time_steps
    steps_done = steps_done + steps;
    liv.steps_done_in_the_memory = int32(double(liv.steps_done_in_the_memory) + steps);
  else
    steps_done = time_steps;
  end
  liv.steps_done = int32(steps_done);

  hashes_now = floor(steps_done / time_steps * double(liv.nhash));
  liv.hashes_now = int32(hashes_now);

  if double(liv.hashes_now) == double(liv.hashes_done)
    tm.live = liv;
    return;
  end

  if isempty(liv.cput_seg) || numel(liv.cput_seg) < 2
    tm = timing.timing_get_time(tm, 'INIT_SEG');
    liv = tm.live;
  end

  seg1 = liv.cput_seg(1, 1);
  etch = timing.timing_string(seg1);
  if isempty(strtrim(etch))
    etch = '--';
  end

  time_report = (seg1 - liv.cput_last_report >= rts) || (liv.steps_done == liv.time_steps);
  if ~time_report
    tm.live = liv;
    return;
  end

  liv.cput_last_report = seg1;

  mems = double(liv.memory_steps);
  if mems > 0
    sdim = double(liv.steps_done_in_the_memory);
    if sdim >= mems
      lts = (seg1 - liv.cput_last_estimate) * (time_steps - steps_done) / max(sdim, 1);
      lts = lts + seg1;
      liv.steps_done_in_the_memory = int32(0);
      liv.cput_last_estimate = seg1;
    else
      lts = seg1 * time_steps / max(steps_done, 1);
    end
  else
    lts = seg1 * time_steps / max(steps_done, 1);
  end

  xtch = timing.timing_string(lts);
  if isempty(strtrim(xtch))
    xtch = '--';
  end

  tm.live = liv;
  tm = live_timing_update(tm, xtch, etch);
end

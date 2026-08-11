% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08
%
% Port of Yambo LIVE_timing.F
%
%   timing.LIVE_timing()                 -> close (flush progress)
%   timing.LIVE_timing(steps)            -> add integer steps (numeric scalar)
%   timing.LIVE_timing(steps, label)     -> add steps and refresh bar with label
%                                          (steps=0 refreshes label only)
%   timing.LIVE_timing(msg, steps)       -> activate with total step count
%   timing.LIVE_timing(msg, steps, depth)-> optional DEPTH for memory_steps
%
% Requires timing.driver first. MATLAB R2008a compatible.

function LIVE_timing(varargin)
  if nargin == 0
    tm = timing.get();
    tm = live_timing_close(tm);
    timing.save2mod(tm);
    return;
  end

  if nargin == 1 && isnumeric(varargin{1})
    tm = timing.get();
    tm = live_timing_add(tm, varargin{1});
    timing.save2mod(tm);
    return;
  end

  % LIVE_timing(steps, label): update bar caption (e.g. relative loss) while advancing.
  if nargin == 2 && isnumeric(varargin{1}) && (ischar(varargin{2}) || isstring(varargin{2}))
    steps = double(varargin{1});
    label = char(string(varargin{2}));
    tm = timing.get();
    if ~tm.live.live_timing_is_on
      return;
    end
    tm.live.timing_name = label;
    hd0 = tm.live.hashes_done;
    if steps > 0 && isfinite(steps) && floor(steps) == steps
      tm = live_timing_add(tm, steps);
    end
    % Refresh when bar did not advance so the new label (e.g. rel loss) still appears.
    if steps <= 0 || tm.live.hashes_done == hd0
      tm = timing.timing_get_time(tm, 'SEG');
      ech = '--';
      if ~isempty(tm.live.cput_seg)
        ech = timing.timing_string(tm.live.cput_seg(1, 1));
        if isempty(strtrim(ech))
          ech = '--';
        end
      end
      tm = live_timing_update(tm, '--', ech, true);
    end
    timing.save2mod(tm);
    return;
  end

  if nargin >= 2
    msg = varargin{1};
    steps = varargin{2};
    depth = [];
    if nargin >= 3
      depth = varargin{3};
    end
    if ~ischar(msg)
      error('timing:LIVE_timing:InvalidMessage', 'Progress label must be char.');
    end
    if ~(isnumeric(steps) && isequal(size(steps), [1, 1]) && isfinite(steps) && steps > 0 && floor(steps) == steps)
      error('timing:LIVE_timing:InvalidSteps', 'steps must be a finite positive integer scalar.');
    end
    tm = timing.get();
    tm = live_timing_activate(tm, msg, int32(steps), depth);
    timing.save2mod(tm);
    return;
  end

  error('timing:LIVE_timing:Syntax', 'Invalid call to LIVE_timing.');
end

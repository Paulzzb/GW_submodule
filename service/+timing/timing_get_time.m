% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/03 ZZ
%
% Wall-clock timing via tic/toc (not cputime). LIVE (E) and section clocks
% report real elapsed time, which is what interactive MATLAB users expect.
% Options (char flags in varargin): 'INIT','INIT_SEG','SEG','INIT_SEC','SEC','FIN'
% Layout: cput_* are 1x2 for compatibility; col1 = elapsed seconds, col2 unused.

function tm = timing_get_time(tm, varargin)
  liv = tm.live;

  if isempty(liv.tic_tot)
    t0 = tic;
    liv.tic_tot = t0;
    liv.tic_sec = t0;
    liv.tic_seg = t0;
    liv.cput_tot = [0, 0];
    liv.cput_sec = [0, 0];
    liv.cput_seg = [0, 0];
  end

  if timing_get_time_has(varargin, 'INIT')
    t0 = tic;
    liv.tic_tot = t0;
    liv.tic_sec = t0;
    liv.tic_seg = t0;
    liv.cput_tot = [0, 0];
    liv.cput_sec = [0, 0];
    liv.cput_seg = [0, 0];
  end

  liv.cput_tot(1, 1) = toc(liv.tic_tot);

  if timing_get_time_has(varargin, 'INIT_SEC')
    liv.tic_sec = tic;
    liv.cput_sec = [0, 0];
  end
  if timing_get_time_has(varargin, 'INIT_SEG')
    liv.tic_seg = tic;
    liv.cput_seg = [0, 0];
  end
  if timing_get_time_has(varargin, 'SEC')
    if isempty(liv.tic_sec)
      liv.tic_sec = tic;
    end
    liv.cput_sec(1, 1) = toc(liv.tic_sec);
  end
  if timing_get_time_has(varargin, 'SEG')
    if isempty(liv.tic_seg)
      liv.tic_seg = tic;
    end
    liv.cput_seg(1, 1) = toc(liv.tic_seg);
  end
  if timing_get_time_has(varargin, 'FIN')
    liv.tic_seg = [];
    liv.tic_sec = [];
    liv.tic_tot = [];
    liv.cput_seg = zeros(0, 0);
    liv.cput_sec = zeros(0, 0);
    liv.cput_tot = zeros(0, 0);
  end

  tm.live = liv;
end

function tf = timing_get_time_has(opts, key)
  tf = false;
  for i = 1:numel(opts)
    v = opts{i};
    if ischar(v) && strcmpi(v, key)
      tf = true;
      return;
    end
  end
end

% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08
%
% Port of Yambo TIMING_get_time.F for serial MATLAB (double "CPU" row).
% Options (char flags in varargin): 'INIT','INIT_SEG','SEG','INIT_SEC','SEC','FIN'
% Layout: cput_* are 1x2, col1 = elapsed, col2 = reference cputime anchor.

function tm = timing_get_time(tm, varargin)
  liv = tm.live;
  cput_now = cputime;

  if isempty(liv.cput_tot) || size(liv.cput_tot, 1) < 1 || size(liv.cput_tot, 2) < 2
    liv.cput_tot = [0, cput_now];
    liv.cput_sec = [0, cput_now];
    liv.cput_seg = [0, cput_now];
  end

  if timing_get_time_has(varargin, 'INIT')
    liv.cput_seg = [0, cput_now];
    liv.cput_sec = [0, cput_now];
    liv.cput_tot = [0, cput_now];
  end

  liv.cput_tot(1, 1) = cput_now - liv.cput_tot(1, 2);

  if timing_get_time_has(varargin, 'INIT_SEC')
    liv.cput_sec(1, :) = [0, cput_now];
  end
  if timing_get_time_has(varargin, 'INIT_SEG')
    liv.cput_seg(1, :) = [0, cput_now];
  end
  if timing_get_time_has(varargin, 'SEC')
    liv.cput_sec(1, 1) = cput_now - liv.cput_sec(1, 2);
  end
  if timing_get_time_has(varargin, 'SEG')
    if isempty(liv.cput_seg) || numel(liv.cput_seg) < 2
      liv.cput_seg = [0, cput_now];
    end
    liv.cput_seg(1, 1) = cput_now - liv.cput_seg(1, 2);
  end
  if timing_get_time_has(varargin, 'FIN')
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

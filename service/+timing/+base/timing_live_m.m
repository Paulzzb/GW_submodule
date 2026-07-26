% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08
%
% Mirrors Yambo services/timing/base/mod_TIMING_live.F
% MATLAB R2008a: char vectors instead of string.

classdef timing_live_m
  properties
    date_time_at_start = zeros(1, 6)
    live_timing_is_on = false

    nhash = int32(0)
    time_steps = int32(0)
    steps_done = int32(0)
    steps_done_in_the_memory = int32(0)
    hashes_now = int32(0)
    hashes_done = int32(0)
    memory_steps = int32(0)

    cput_seg = zeros(0, 0)
    cput_sec = zeros(0, 0)
    cput_tot = zeros(0, 0)

    cput_last_report = 0
    cput_last_estimate = 0
    timing_name = ''

    USER_wall_time_string = ' '
    USER_wall_time = [0, 0, 0]

    log_line_to_dump = false
    log_line = ' '

    % Minimum segment CPU time (seconds) between LIVE bar prints when the hash
    % advances; 0 => print on every hash change (closer to isdf_coarse_validate_energies).
    % Yambo LIVE_timing_add uses rts = 5.
    live_report_min_seconds = 5

    % Empty => timing.report writes to command window; otherwise append path (char).
    report_logfile = ''
  end

  methods
    function obj = timing_live_m()
      % Defaults aligned with TIMING_defaults.F / mod_TIMING_live.F
      obj.nhash = int32(40);
    end
  end
end

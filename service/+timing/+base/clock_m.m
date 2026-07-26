% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08
%
% Mirrors Yambo services/timing/base/mod_TIMING.F TYPE(clock).
% Written for MATLAB R2008a-style OOP (char, no property validation).

classdef clock_m
  properties
    name = ''
    cpu_id = int32(-1)
    call_number = int32(0)
    start = 0
    stop = 0
    total_time = 0
    running = false
    indx = int32(0)
    alloc = false
  end

  methods
    function obj = clock_m()
    end
  end
end

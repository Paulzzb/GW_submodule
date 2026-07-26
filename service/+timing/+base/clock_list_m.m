% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08
%
% Mirrors Yambo TYPE(clock_list) + CLOCK_list_allocate / CLOCK_allocate.
% MATLAB R2008a: no validateattributes / string / .empty on user classes.

classdef clock_list_m
  properties
    name = ''
    nclock = int32(0)
    nclock_max = int32(0)
    alloc = false
    clocks = []
  end

  methods
    function obj = clock_list_m()
    end
  end

  methods (Static)
    function list = preallocate(list_name, nclock_max_)
      if ~(isnumeric(nclock_max_) && isequal(size(nclock_max_), [1, 1]) && isfinite(nclock_max_) && nclock_max_ > 0 && floor(nclock_max_) == nclock_max_)
        error('clock_list_m:preallocate:InvalidNclockMax', 'nclock_max_ must be a finite positive integer scalar.');
      end

      nm = strtrim(char(list_name));
      if isempty(nm)
        error('clock_list_m:preallocate:InvalidName', 'List name must be non-empty.');
      end

      list = timing.base.clock_list_m;
      list.name = nm;
      list.nclock_max = int32(nclock_max_);
      list.nclock = int32(0);
      list.clocks = repmat(timing.base.clock_m(), double(nclock_max_), 1);
      list.alloc = true;
    end
  end

  methods
    function list = allocate_next_clock(list, clock_name)
      if ~list.alloc
        error('clock_list_m:allocate_next_clock:NotAllocated', 'Clock list is not allocated.');
      end

      nm = strtrim(char(clock_name));
      if isempty(nm)
        error('clock_list_m:allocate_next_clock:InvalidName', 'Clock name must be non-empty.');
      end

      idx = double(list.nclock) + 1;
      if idx > double(list.nclock_max)
        error('clock_list_m:allocate_next_clock:TooManyClocks', 'Too many clocks in list ''%s''.', list.name);
      end

      c = list.clocks(idx);
      if c.alloc
        error('clock_list_m:allocate_next_clock:AlreadyAllocated', 'Clock slot %d already allocated.', idx);
      end

      c.name = nm;
      c.cpu_id = int32(-1);
      c.call_number = int32(0);
      c.start = 0;
      c.stop = 0;
      c.total_time = 0;
      c.running = false;
      c.indx = int32(idx);
      c.alloc = true;

      list.clocks(idx) = c;
      list.nclock = int32(idx);
    end
  end
end

% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08
%
% Mirrors Yambo services/timing/base/mod_TIMING.F module state
% (internal_list, global_list, verb, overview widths, alloc) plus
% embedded TIMING_live / TIMING_logo payloads as MATLAB objects.
%
% Style targets MATLAB R2008a: char (not string), no property validation
% attributes, no validateattributes; see also timing.driver (addParamValue).

classdef timing_m
  properties (Constant)
    NCLOCKX = int32(200)
  end

  properties
    internal_list = timing.base.clock_list_m()
    global_list = timing.base.clock_list_m()

    TIMING_verb = int32(0)

    max_name_length = int32(0)
    max_calls_length = int32(0)

    alloc = false

    live = timing.base.timing_live_m()
    logo = timing.base.timing_logo_m()
  end

  methods (Static)
    function obj = create_default(nclock_max_global)
      if nargin < 1 || isempty(nclock_max_global)
        nmax = timing.base.timing_m.NCLOCKX;
      else
        if ~(isnumeric(nclock_max_global) && isequal(size(nclock_max_global), [1, 1]) && isfinite(nclock_max_global) && nclock_max_global > 0 && floor(nclock_max_global) == nclock_max_global)
          error('timing_m:create_default:InvalidNclockMax', 'nclock_max_global must be a finite positive integer scalar.');
        end
        nmax = int32(nclock_max_global);
      end

      obj = timing.base.timing_m;
      obj.global_list = timing.base.clock_list_m.preallocate('global', nmax);
      obj.internal_list = timing.base.clock_list_m.preallocate('internal', int32(1));
      obj.internal_list = obj.internal_list.allocate_next_clock('internal');
      obj.live = timing.base.timing_live_m();
      obj.logo = timing.base.timing_logo_m();
      obj = obj.apply_logo_defaults();
      obj.alloc = true;
    end
  end

  methods (Access = private)
    function obj = apply_logo_defaults(obj)
      % TIMING_defaults.F: ID_logo / ID_logo_stderr = -1
      lg = obj.logo;
      lg.ID_logo = int32(-1);
      lg.ID_logo_stderr = int32(-1);
      lg.n_logo_lines = int32(0);
      obj.logo = lg;
    end
  end

  methods
    function obj = timing_m()
    end
  end
end

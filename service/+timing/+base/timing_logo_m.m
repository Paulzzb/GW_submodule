% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08
%
% Mirrors Yambo services/timing/base/mod_TIMING_logo.F
% MATLAB R2008a: cell column of char rows (Fortran character(70) lines).

classdef timing_logo_m
  properties (Constant)
    % Fortran max_n_logo_lines=100 (do not read this Constant from the ctor:
    % timing_m default property init cannot resolve it in some MATLAB releases.)
    MAX_N_LOGO_LINES = int32(100)
  end

  properties
    n_logo_lines = int32(0)
    ID_logo = int32(-1)
    ID_logo_stderr = int32(-1)
    logo_line = {}
  end

  methods
    function obj = timing_logo_m()
      n = 100;
      obj.logo_line = repmat({''}, n, 1);
    end
  end
end

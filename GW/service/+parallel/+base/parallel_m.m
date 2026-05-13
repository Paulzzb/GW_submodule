% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/06

classdef parallel_m
  properties
    enabled(1, 1) logical = false
    requested(1, 1) logical = true
    toolbox_available(1, 1) logical = false
    workers(1, 1) {mustBeInteger} = int32(0)
    mode(1, 1) string = "serial"
  end
end

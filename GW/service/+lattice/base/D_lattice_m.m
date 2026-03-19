% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/18 ZZ

classdef D_lattice_m
  properties
    a1a2a3(3, 3) single = single(zeros(3, 3))  % lattice vectors (columns)
    DL_vol(1, 1) single = single(0)  % unit cell volume
    atom_pos(:, :, 3) single = single(zeros(0, 0, 3))  % atomic positions (Cartesian)
  end
  
  methods
    function obj = D_lattice_m()
      % Constructor (currently does nothing)
      obj.a1a2a3 = single(zeros(3, 3));
      obj.DL_vol = single(0);
      obj.atom_pos = single(zeros(0, 0, 3));
    end
  end
end

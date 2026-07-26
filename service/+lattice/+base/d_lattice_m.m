% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/18 ZZ

classdef d_lattice_m
  properties
    a1a2a3(3, 3) double = double(zeros(3, 3))  % lattice vectors (columns)
    DL_vol(1, 1) double = double(0)  % unit cell volume
    atom_pos(:, :, 3) double = double(zeros(0, 0, 3))  % atomic positions (Cartesian)
    atom_symbol(:, 1) cell = cell(0, 1)  % per-atom element names (e.g. {'Si','O'})
  end
  
  methods
    function obj = d_lattice_m()
      % Constructor (currently does nothing)
      obj.a1a2a3 = double(zeros(3, 3));
      obj.DL_vol = double(0);
      obj.atom_pos = double(zeros(0, 0, 3));
      obj.atom_symbol = cell(0, 1);
    end
  end
end

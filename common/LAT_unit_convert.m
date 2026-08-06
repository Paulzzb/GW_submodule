% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/19 ZZ

function v_out = LAT_unit_convert(v1v2v3, v_in, mode)
% Coordinate conversion used in Yambo
%
% v_out = LAT_unit_convert(b_in, v_in, mode)
%
% INPUT
%   b_in : 3x3 lattice matrix
%   v_in : 1x3 vector
%   mode : string
%
% OUTPUT
%   v_out : converted vector
%
% Required global variables (same logic as Yambo)
%   alat : lattice constants [1x3]
%   a    : direct lattice vectors (3x3)
%   b    : reciprocal lattice vectors (3x3)

% global alat a b

% ---- choose lattice ----
a_here = v1v2v3;
alat = a_here(1:3,1);

% ---- scaling factors ----
n = ones(1,3);

if contains(mode,'r')
  n = alat;
end

if contains(mode,'k')
  n = 2*pi ./ alat;
end

% ---- convert iku -> cc ----
u = v_in;

if contains(mode,'i2')
  u = v_in .* n;
end

% ---- transformation matrix ----
mat = eye(3);

% a2c / a2i
if contains(mode,'a2c') || contains(mode,'a2i')
  mat = a_here.';
end

% c2a / i2a
if contains(mode,'c2a') || contains(mode,'i2a')
  mat = inv(a_here.');
end

% ---- apply transform ----
v_out = (mat * u.').';

% ---- convert to iku ----
if contains(mode,'2i')
  v_out = v_out ./ n;
end

end
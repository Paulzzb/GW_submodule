% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/18 ZZ

classdef FFT_m
  properties
    fftgrid(1, 3) {mustBeInteger} = int32(zeros(1, 3))
    nr(1, 1) {mustBeInteger} = int32(0)
    % ng(1, 1) {mustBeInteger} = int32(0)
    Rgrid_RLU(:, 3) {mustBeInteger} = int32(zeros(0, 3))
    R_rot(:, :) {mustBeInteger} = int32(zeros(0, 0))
    R_rot_inv(:, 1) {mustBeInteger} = int32(zeros(0, 1))
    nGo(1, 1) {mustBeInteger} = int32(0)
    G_table(:, :) {mustBeInteger} = int32(zeros(0, 0))
  end

  methods
    function obj = FFT_m(fftgrid, nsym)
      if ~(isequal(size(fftgrid), [1, 3]) || isequal(size(fftgrid), [3, 1]))
        error('FFT_m:InvalidFFTGridShape', 'fftgrid must be a 1x3 or 3x1 integer vector.');
      end

      validateattributes(fftgrid, {'numeric'}, {'real', 'finite', 'integer', 'positive'}, mfilename, 'fftgrid');
      validateattributes(nsym, {'numeric'}, {'real', 'finite', 'scalar', 'integer', 'positive'}, mfilename, 'nsym');

      fftgrid = int32( reshape(fftgrid, 1, 3) );
      nsym = int32(nsym);
      ngrid_pts = prod(double(fftgrid));

      obj.fftgrid = fftgrid;
      obj.nr = int32(ngrid_pts);
      obj.Rgrid_RLU = int32( zeros(ngrid_pts, 3) );
      obj.R_rot = int32( zeros(ngrid_pts, nsym) );
      % obj.G_rot = int32( zeros(ng, nsym) );
      % obj.G_table = int32( zeros(ngrid_pts, nsym) );
    end
  end
end
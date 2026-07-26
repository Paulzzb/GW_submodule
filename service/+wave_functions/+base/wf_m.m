% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/24

classdef wf_m
  properties
    assigned(1, 1) logical = false
    allocated(1, 1) logical = false
    % what(1, 1) string = ""

    % Dynamic wave-function metadata.
    ng(1, 1) {mustBeInteger} = int32(0)
    %
    space (1, 1) string = "r"
    nc(1, 1) {mustBeInteger} = int32(0)
    nb(1, 1) {mustBeInteger} = int32(0)
    nk(1, 1) {mustBeInteger} = int32(0)
    nspin(1, 1) {mustBeInteger} = int32(0)
    sp_pol(1, 2) {mustBeInteger} = int32(zeros(1, 2))
    n_spinor(1, 1) {mustBeInteger} = int32(1)
    c(:, :, :, :) double = complex(zeros(0, 0, 0, 0, 'double'))
  end

  methods
    function obj = wf_m(nc, nb, nk, nspin)
    
      if nargin == 0
        return
      end

      if nargin ~= 4
        error('wf_m:InvalidConstructor', ...
          'Expected constructor arguments (nc, nspin, nstates, index_shape) or no arguments.');
      end

      validateattributes(nc, {'numeric'}, {'scalar', 'integer', 'nonnegative'}, mfilename, 'nc');
      validateattributes(nb, {'numeric'}, {'scalar', 'integer', 'positive'}, mfilename, 'nb');
      validateattributes(nk, {'numeric'}, {'scalar', 'integer', 'nonnegative'}, mfilename, 'nk');
      validateattributes(nspin, {'numeric'}, {'scalar', 'integer', 'positive'}, mfilename, 'nspin');

      obj.nc = int32(nc);
      obj.nb = int32(nb);
      obj.nk = int32(nk);
      obj.nspin = int32(nspin);
      obj.c = complex(zeros(nc, nb, nk, nspin, 'double'));
      obj.n_spinor = int32(1);
      obj.allocated = true;
    end
  end
end
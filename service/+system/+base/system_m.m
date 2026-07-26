% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/26

classdef system_m
  properties
    assigned(1, 1) logical = false
    allocated(1, 1) logical = false

    nb(1, 1) {mustBeInteger} = int32(0)
    nk(1, 1) {mustBeInteger} = int32(0)
    nspin(1, 1) {mustBeInteger} = int32(0)

    Eo(:, :, :) double = zeros(0, 0, 0)
    Vxc(:, :, :) double = zeros(0, 0, 0)
    Eqp(:, :, :) double = zeros(0, 0, 0)
    qptype (1, 1) string = "HF"
    f(:, :, :) double = zeros(0, 0, 0)
    first_index_in_degeneracy(:, 1) cell = cell(0, 1)
    num_index_in_degeneracy(:, 1) cell = cell(0, 1)
    degeneracy_indices_len(:, 1) int32 = zeros(0, 1, 'int32')
  end

  methods
    function obj = system_m(nb, nk, nspin)
      if nargin == 0
        return
      end

      if nargin ~= 3
        error('system_m:InvalidConstructor', ...
          'Expected constructor arguments (nb, nk, nspin) or no arguments.');
      end

      validateattributes(nb, {'numeric'}, {'scalar', 'integer', 'nonnegative'}, mfilename, 'nb');
      validateattributes(nk, {'numeric'}, {'scalar', 'integer', 'nonnegative'}, mfilename, 'nk');
      validateattributes(nspin, {'numeric'}, {'scalar', 'integer', 'positive'}, mfilename, 'nspin');

      obj.nb = int32(nb);
      obj.nk = int32(nk);
      obj.nspin = int32(nspin);
      obj.Eo = zeros(nb, nk, nspin);
      obj.Vxc = zeros(nb, nk, nspin);
      obj.Eqp = zeros(nb, nk, nspin);
      obj.qptype = "HF";
      obj.f = zeros(nb, nk, nspin);
      obj.first_index_in_degeneracy = cell(nk, 1);
      obj.num_index_in_degeneracy = cell(nk, 1);
      obj.degeneracy_indices_len = zeros(nk, 1, 'int32');
      obj.allocated = true;
    end
  end
end
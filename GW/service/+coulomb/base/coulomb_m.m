% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/18 ZZ

classdef coulomb_m
  properties
    % Placeholder for properties related to Coulomb matrix
    % For example, you might have:
    % coulomb_matrix(:, :) {mustBeNumeric} = zeros(0, 0)
    trunc_method(1, 1) {mustBeInteger} = int32(0)
    % 0: no truncation, see ./doc for detail
    trunc_param(1, 1) {mustBeNumeric} = 0
    %
    coulomb_ng(1, 1) {mustBeInteger} = int32(0)
    bare_qpg(:, :) = zeros(0, 0)
    %
    vcoul(:, :) = zeros(0, 0)
    vcoul0(1, 1) = zeros(0, 0) % at q=0, G=0
    %
  end

  methods
    function obj = coulomb_m(ng, nqibz)
      % Constructor for coulomb_m class
      % Initialize properties as needed
      validateattributes(ng, {'numeric'}, {'integer', 'positive'}, mfilename, 'ng');
      validateattributes(nqibz, {'numeric'}, {'integer', 'positive'}, mfilename, 'nqibz');
      obj.coulomb_ng = int32(ng);
      obj.bare_qpg = zeros(ng, nqibz);
      obj.vcoul = zeros(ng, nqibz);
    end
  end
end

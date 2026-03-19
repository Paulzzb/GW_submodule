% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/18 ZZ

classdef bz_samp_m
  properties
    %------------------------------------------------------------
    % Coarse grid (IBZ + full BZ)
    %------------------------------------------------------------
    description(1, :) char = ''   % e.g. 'k, q'
    %
    nibz(1, 1)  {mustBeInteger, mustBeNonnegative} = int32(0)  % # of irreducible k-points
    nbz(1, 1)   {mustBeInteger, mustBeNonnegative} = int32(0)  % # of full-BZ k-points

    % nstar(nibz)  : number of BZ k-points in the star of each IBZ k-point
    nstar(:, 1) {mustBeInteger} = int32(zeros(0, 1))

    % star(nibz, nsym) : ik_ibz, i_star -> isym that sends ik_ibz to the i_star-th BZ point
    star(:, :)  {mustBeInteger} = int32(zeros(0, 0))

    % sstar(nbz, 2) : ik_bz -> [ik_ibz, isym]  (inverse of star)
    sstar(:, 2) {mustBeInteger} = int32(zeros(0, 2))

    % s_table(nibz, nsym) : ik_ibz, isym -> ik_bz  (reverse of sstar(:,2))
    s_table(:, :) {mustBeInteger} = int32(zeros(0, 0))

    % k_table(nibz, nsym) : ik_ibz, isym -> ik_bz  (reverse of sstar(:,1))
    k_table(:, :) {mustBeInteger} = int32(zeros(0, 0))


    %------------------------------------------------------------
    % Extended fields used by the GW codebase
    %------------------------------------------------------------
    b1b2b3(3, 3)   = single(zeros(3, 3))  % Reciprocal lattice matrix (rows = b1,b2,b3)

    % kpt(nibz, 3)   : RLU k-points  (alias for pt, Cartesian)
    kpt_RLU(:, 3)     = single(zeros(0, 3))
    kpt_Cart(:, 3)     = single(zeros(0, 3))
    % weights(nibz) : integration weights for IBZ k-points
    weights(:, 1) = single(zeros(0, 1))

    % kptbz(nbz, 3)  : full-BZ k-points (alias for ptbz, Cartesian)
    kptbz_RLU(:, 3)   = single(zeros(0, 3))
    kptbz_Cart(:, 3)   = single(zeros(0, 3))

    bz2ibz(:, 1)        {mustBeInteger} = int32(zeros(0, 1))  % ik_bz -> ik_ibz
    bz2rot(:, 1)        {mustBeInteger} = int32(zeros(0, 1))  % ik_bz -> isym 

    %------------------------------------------------------------
    % Fine grids (reserved; not yet implemented)
    %------------------------------------------------------------
    % FGbare, FGibz, FGbz : bz_fine_grid objects (to be added when needed)
  end

  methods
    function obj = bz_samp(nibz, nsym)
      % bz_samp(nibz, nsym)
      %   Allocate all integer index arrays given the grid sizes.
      %   Floating-point arrays (kpt, kptbz, weights, ...) are left as
      %   empty placeholders and should be filled by the driver.
      if nargin == 0
        return;  % default construction with empty arrays
      end

      nbz = 0;

      validateattributes(nibz, {'numeric'}, {'scalar', 'integer', 'nonnegative'}, mfilename, 'nibz');
      validateattributes(nsym, {'numeric'}, {'scalar', 'integer', 'positive'},    mfilename, 'nsym');

      obj.nibz  = int32(nibz);
      obj.nbz   = int32(nbz);
      obj.nstar = int32(zeros(nibz, 1));
      obj.star  = int32(zeros(nibz, nsym));
      obj.sstar = int32(zeros(nbz,  2));
      obj.s_table = int32(zeros(nibz, nsym));
      obj.k_table = int32(zeros(nibz, nsym));
      obj.kpt_RLU   = single(zeros(nibz, 3));
      obj.kpt_Cart = single(zeros(nibz, 3));
      obj.kptbz_RLU   = single(zeros(nbz, 3));
      obj.kptbz_Cart = single(zeros(nbz, 3));
      obj.kpt_weights = single(zeros(nibz, 1));
      obj.bz2ibz = int32(zeros(nbz, 1));
      obj.bz2rot = int32(zeros(nbz, 1));
      obj.qindx_X = int32(zeros(nibz, nbz, 2));
      obj.qindx_S = int32(zeros(nibz, nbz, 2));
      obj.qindx_C = int32(zeros(nbz, nbz, 2));
      obj.qindx_B = int32(zeros(nbz, nbz, 2));
    end % constructor
  end % methods
end % classdef bz_samp

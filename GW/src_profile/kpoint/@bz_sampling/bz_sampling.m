classdef bz_sampling
  % bz_sampling: Brillouin zone sampling information and symmetry mapping
  % Converted from Fortran type(bz_samp)

  properties
    % --- Coarse Grid ---
    nibz (1,1) {mustBeInteger} = 0           % # of points in irreducible BZ
    nbz  (1,1) {mustBeInteger} = 0           % # of points in full BZ

    nstar (:,1)   % [nibz x 1] number of points in symmetry star
    star  (:,:)   % [nibz x max_star] maps ik,ik* -> is
    sstar (:,:)   % [nbz x 2] maps ik_bz -> [ik_ibz, is]
    s_table (:,:) % [nibz x max_star] reverse of sstar(2,:)
    k_table (:,:) % [nibz x max_star] reverse of sstar(1,:)

    kpt    (:,:)   % [3 x nibz] fractional coordinates of IBZ points
    kptbz  (:,:)   % [3 x nbz]  fractional coordinates of full BZ points
    weights (:,1) % [nibz x 1] weights of IBZ points
    bmatrix (3, 3) % [3 x 3] reciprocal lattice vectors

    description (1,:) char = ""  % e.g. 'MP grid'
    units       (1,1) char = 'r' % 'r': reciprocal; 'c': cartesian

    % --- Fine Grids (optional substructure) ---
    FGbare   % can be another class or struct
    FGibz
    FGbz

    % --- Optional for Wannier interpolation ---
    weights_ipol (:,1)  % mapping weights (e.g. for wannier interpolation)
  end

  methods
    function obj = bz_sampling(nibz, kibz, bmatrix)
      obj.nibz = nibz;
      obj.kpt = kibz;
      obj.bmatrix = bmatrix;
    end

    function show(obj)
      fprintf('bz_sampling: %d IBZ points, %d BZ points\n', obj.nibz, obj.nbz);
      fprintf('Description: %s\n', obj.description);
    end
  end
end

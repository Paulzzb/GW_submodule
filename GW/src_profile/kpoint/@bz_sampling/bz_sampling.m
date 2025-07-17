classdef bz_sampling
  % bz_sampling: Brillouin zone sampling information and symmetry mapping
  % Converted from Fortran type(bz_samp)

  properties
    % --- Coarse Grid ---
    nibz (1,1) {mustBeInteger} = 0           % # of points in irreducible BZ
    nbz  (1,1) {mustBeInteger} = 0           % # of points in full BZ

    kpt          (:,:)   % [nbz x 3] fractional coordinates of IBZ points
    kptbz        (:,:)   % [nbz x 3]  fractional coordinates of full BZ points
    kptweights   (:,1)  % [nibz x 1] weights of IBZ points
    bmatrix (3, 3) % [3 x 3] reciprocal lattice vectors

    % Rotation{rotation} * kptbz(ind_kbz) = kpt(ind_rotation) + ind_G0(ind_G0)
    kbz2kibz_ind_rotation (:, 1) {mustBeInteger} % [nbz x 1]
    kbz2kibz_ind_kbz      (:, 1) {mustBeInteger} % [nbz x 1] 
    kbz2kibz_ind_Go       (:, 1) {mustBeInteger} % [nbz x 1] 

    % k-q --> k' mapping
    nstar        (:,1)   % [nibz x 1] number of points in symmetry star
    star         (:,:)   % [nibz x max_star] maps ik,ik* -> is
    sstar        (:,:)   % [nbz x 2] maps ik_bz -> [ik_ibz, is]
    s_table      (:,:) % [nibz x max_star] reverse of sstar(2,:)
    k_table      (:,:) % [nibz x max_star] reverse of sstar(1,:)
    
    % qindx_S(ik,iqbz,1)=okbz
    % qindx_S(ik,iqbz,2)=iGo
    qindx_S      (:,:,:) {mustBeInteger} % [nibz x nbz x 2] 
    % qindx_C(ikbz,iqbz,1)=okbz
    % qindx_C(ikbz,iqbz,2)=iGo
    qindx_C      (:,:,:) {mustBeInteger} % [nbz x nbz x 2] 
    % qindx_X(iq,ikbz,1)=okbz
    % qindx_X(iq,ikbz,2)=iGo
    qindx_X      (:,:,:) {mustBeInteger} % [nibz x nbz x 2] 
    % qindx_B(okbz,ikbz,1)=iqbz
    % qindx_B(okbz,ikbz,2)=iGo
    qindx_B      (:,:,:) {mustBeInteger} % [nbz x nbz x 2] 


    % description (1,:) char = ""  % e.g. 'MP grid'
    % units       (1,1) char = 'r' % 'r': reciprocal; 'c': cartesian

    

    % Possible G0set and corresponding index
    

    % Possible G0set and corresponding index
    nGo          (1,1) {mustBeInteger} = 1 
    iGolist      (:,1) {mustBeInteger} = 0 % [nG0 x 3] G0 vectors

  end

  methods
    function obj = bz_sampling(nibz, kibz, bmatrix, weights)
      obj.nibz = nibz;
      obj.kpt = kibz;
      obj.bmatrix = bmatrix;
      obj.kptweights = weights;
    end

    function show(obj)
      fprintf('bz_sampling: %d IBZ points, %d BZ points\n', obj.nibz, obj.nbz);
      fprintf('Description: %s\n', obj.description);
    end
  end
end

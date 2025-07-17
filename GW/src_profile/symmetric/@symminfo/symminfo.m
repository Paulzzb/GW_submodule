classdef symminfo
  % symminfo: Encapsulates symmetry operations for BZ sampling
  
  properties
    ntran   (1,1) {mustBeInteger} = 0    % # symmetry ops preserving q
    ntranq  (1,1) {mustBeInteger} = 0    % unused (reserved)
    mtrx    cell                         % {nrot x 1} rotation matrices
    nrot    (1,1) {mustBeInteger} = 0    % total # symmetry operations
    % SURE??
    indsub  (:,1) double  = 0            % sub-operation index
    kgzero  (:,3) double  = 0            % G vectors left invariant
    % 
    % ng      (1,1) {mustBeInteger} = 0    % number of G vectors, = max_ng from all gvec
    % R_is G_ig = G_{g_rot(ig,is)} 
    % g_rot   (:, :) {mustBeInteger}       % [ng x nrot]  
  end

  methods
    function obj = symminfo(ntran, ntranq, mtrx, nrot, indsub, kgzero)
      if nargin > 0
        obj.ntran  = ntran;
        obj.ntranq = ntranq;
        obj.mtrx   = mtrx;
        obj.nrot   = nrot;
        obj.indsub = indsub;
        obj.kgzero = kgzero;
      end
    end

    function show(obj)
      fprintf('Symmetry info: %d total rotations, %d invariant under q\n', ...
              obj.nrot, obj.ntran);
    end

    function R = get_rotation(obj, i)
      % Return the i-th rotation matrix
      R = obj.mtrx{i};
    end

    function Gp = apply_symmetry(obj, i, G)
      % Apply i-th symmetry op to G-vector: Gp = R * G
      R = obj.get_rotation(i);
      Gp = (R * G')';  % Apply as row vector
    end
  end
end

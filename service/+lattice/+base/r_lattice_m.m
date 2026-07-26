classdef r_lattice_m
  properties
    b1b2b3(3, 3) double = double(zeros(3, 3)) % columns are b1, b2, b3
    RL_vol(1, 1) double = double(0) % volume of the reciprocal lattice
    d3k_factor(1, 1) double = double(0) % prefactor for k-point integration: (2pi)^3 / RL_vol
    d3q_factor(1, 1) double = double(0) % prefactor for q-point integration: (2pi)^3 / RL_vol
    
    k % bz_samp_m object for k-points
    q % bz_samp_m object for q-points

    % qindx_S(nibz, nbz, 2) : [ik_kp_bz, iG0] for x=sigma matrix element
    qindx_S = int32(zeros(0, 0, 2))
    qindx_X = int32(zeros(0, 0, 2))
    qindx_C = int32(zeros(0, 0, 2))
    qindx_B = int32(zeros(0, 0, 2))
    qindx_S_max_Go(1, 1) int32 = int32(0) % maximum G-shell index for qindx_S

    % RL_lattice
    ng = int32(0) % number of G-vectors within the cutoff
    Ggrid_RLU(:, 3) int32 = int32( zeros(0, 3) ) % G-vectors in RLU
    Ggrid_Cart(:, 3) double = double(zeros(0, 3)) % G-vectors in Cartesian coordinates
    idxnz(:, 1) int32 = int32(zeros(0, 1)) % indices of G-vectors within the cutoff in the full FFT grid
    %
    n_g_shell(1, 1) int32 = int32(0) % number of G-shells
    num_index_in_each_shell(:, 1) int32 = int32(zeros(0, 1)) % number of G-vectors in each shell
    first_index_in_each_shell(:, 1) int32 = int32(zeros(0, 1)) % first index in each G-shell
    % G_vec(ig, irot) = S(irot) * G_vec(ig)
    G_rot(:, :) int32 = int32(zeros(0, 0))
    tol = 1e-3 % tolerance for matching k-points in qindx
  end
  methods
    function obj = r_lattice_m()
    end
  end

end

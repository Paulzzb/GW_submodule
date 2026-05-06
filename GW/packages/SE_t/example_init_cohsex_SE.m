% EXAMPLE_INIT_COHSEX_SE  Build one SE_t object for static CoHSEX (serial)
%
% Run from this folder or addpath('.../matlab_SE_t').

clc
se = SE_t();
se.touch_assigned();
se.active = true;
se.kind = 'Self-Energy';
se.approx = 'cohsex';
se.XC_id = int32(0);
se.NW = int32(1);
se.NG = int32(500);          % example G-shell count
se.NB = int32([1, 10]);      % valence band range for GW_driver inner loop
se.prefactor = single(1);
se.from_DB = false;
se.desc = struct();         % no DESC_t port

% Frequency mesh (static)
se.W = yambo_w_samp_static();

% Basis: NS_todo transitions (n, n', k) — columns [ib, ibp, ik_ibz]
table_nk = int32([
    1 1 1
    2 2 1
    3 3 1
    ]);
se.basis = yambo_basis_minimal(size(table_nk, 1), table_nk, 'NK', int32(1));

% Parallel schemes: everyone owns everything
n_qp = double(se.basis.NS_todo);
n_b = double(se.NB(2) - se.NB(1) + 1);
n_q = int32(8);              % example q-grid size; GW uses q%nbz

se.PAR_qp = yambo_serial_par_scheme(n_qp, n_qp);
se.PAR_b = yambo_serial_par_scheme(n_b, n_b);
se.PAR_RL = yambo_serial_par_scheme(1, 1);
se.PAR_q = yambo_serial_par_scheme(double(n_q), double(n_q));
se.PAR_structure = int32(1);

% Self-energy accumulator SF(state, freq) — COHSEX uses NW=1
se.alloc_SF(double(se.basis.NS_todo), double(se.NW));

disp('SE_t (MATLAB) summary:')
disp(se)

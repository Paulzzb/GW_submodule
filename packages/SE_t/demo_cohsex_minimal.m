%% DEMO_COHSEX_MINIMAL  Toy COHSEX (static) walk-through without SCATTER_*
%
% Mirrors the numerical core of:
%   GW_COHSEX_residuals  ->  M <- M .* gamp
%   GW_COHSEX_accumulate ->  pre_factor = Re/contract(rho, M*conj(rho));
%                            SEX: -= (4/spin_occ)*pi*pre_factor*f
%                            COH: += 2*pi*pre_factor
%
% Run: addpath(fileparts(mfilename('fullpath'))); demo_cohsex_minimal

%% -------------------------------------------------------------------------
%% PORT-NOTE [SCATTER_* — not called in this demo]
%
% In Fortran (packages/GW/GW_driver.F) the full pipeline uses:
%
%   SCATTER_Gamp_gpu(iq_or_1, isc, 'c')
%       % fills isc%gamp(G,G') from Coulomb / (eps^{-1}) head (see SCATTER_Gamp.F)
%
%   SCATTER_GW_kinematics(iq, i_bas, q, k, SE(ID)%basis, isc, iscp)
%       % sets isc%is, isc%os, isc%qs (band / k / q indices only)
%
%   SCATTER_Bamp_gpu(isc)   % or (ID_wf, ID_fft, isc) in other drivers
%       % builds isc%rhotw(G) = <ib,k|e^{iG r}|ob,k-q> via WF + FFT
%
% Here we inject toy gamp and rho instead.
%% -------------------------------------------------------------------------

addpath(fileparts(mfilename('fullpath')));

%% --- Parameters (toy) -----------------------------------------------------
ng        = 32;                 % G-vector count (toy)
spin_occ  = 2;                  % 2 = spin-saturated / closed shell style (Fortran spin_occ)
n_ibz     = 1;                  % one k-point in toy world
n_qbz     = 4;                  % toy q-shell count

Se = SE_t();
Se.touch_assigned();
Se.active   = true;
Se.kind     = 'Self-Energy';
Se.approx   = 'cohsex';
Se.NW       = int32(1);
Se.NG       = int32(ng);
Se.NB       = int32([1, 4]);    % inner band loop (occupied-like)
Se.prefactor = double(1);

table_nk = int32([
    1 1 1
    2 2 1
    ]);
Se.basis = yambo_basis_minimal(size(table_nk, 1), table_nk, 'NK', int32(n_ibz));
Se.W = yambo_w_samp_static();

n_states = double(Se.basis.NS_todo);
n_b      = double(Se.NB(2) - Se.NB(1) + 1);
Se.PAR_qp = yambo_serial_par_scheme(n_states, n_states);
Se.PAR_b  = yambo_serial_par_scheme(n_b, n_b);
Se.PAR_q  = yambo_serial_par_scheme(n_qbz, n_qbz);
Se.PAR_RL = yambo_serial_par_scheme(1, 1);

Se.alloc_SF(n_states, 1);

%% --- Toy "electrons": occupation f(ib, ik, isp) ---------------------------
% Algebra check (same as GW_COHSEX_accumulate.F):
%   SEX:  -= (4/spin_occ)*pi*pre*f
%   COH:  += 2*pi*pre
% If spin_occ==2 and f==1, the two cancel: net += 2*pi*pre*(1 - f) = 0.
% So use partial occupation below so SF is visibly non-zero in this toy demo.
E_f = zeros(max(Se.NB(2), 8), n_ibz, 1, 'double');
for ib = Se.NB(1):Se.NB(2)
    E_f(ib, 1, 1) = double(0.65);   % 0 < f < 1  =>  non-zero (2*pi*pre*(1 - 2*f/spin_occ*...))
end

%% --- Main loop skeleton (matches GW_driver shape, SCATTER-free) -----------
for i_count = 1:double(Se.basis.NS_todo)
    if ~Se.PAR_qp.IND.element_1D(i_count)
        continue;
    end
    i_bas = double(Se.basis.states_todo(i_count));

    DP_Sc = zeros(2, 1);   % DP_Sc(1)=SEX track, DP_Sc(2)=COH track (Fortran)

    % ---------------------------------------------------------------------
    % PORT-NOTE [SCATTER_Gamp_gpu]: first call at iq reference
    %   call SCATTER_Gamp_gpu(1, isc, 'c')
    % isc%iqref = 1  ...  when q changes:
    %   call SCATTER_Gamp_gpu(isc%qs(2), isc, 'c')
    % ---------------------------------------------------------------------
    % Toy: same random-ish gamp for every iq (replace with SCATTER output)
    rng(1 + i_bas, 'twister');
    gamp = complex(randn(ng, ng, 'double'), randn(ng, ng, 'double')) * double(0.02);
    gamp = gamp + gamp.';  % symmetrize a bit (not physical, demo only)

    for iq = 1:n_qbz
        % -----------------------------------------------------------------
        % PORT-NOTE [SCATTER_GW_kinematics]
        %   call SCATTER_GW_kinematics(iq, i_bas, q, k, SE(ID)%basis, isc, iscp)
        % sets isc%os(1) later in band loop; isc%os(2)=ik, isc%os(4)=isp etc.
        % -----------------------------------------------------------------
        if ~Se.PAR_q.IND.element_1D(iq)
            continue;
        end

        % Toy static polarizability block chi0(G,G') — from DB / X in real code
        M_chi = complex(randn(ng, ng, 'double'), randn(ng, ng, 'double')) * double(0.01);

        % === GW_COHSEX_residuals (Fortran): M <- M .* gamp =================
        M_eff = cohsex_residuals_toy(M_chi, gamp);

        for ib = double(Se.NB(1)):double(Se.NB(2))
            if ~Se.PAR_b.IND.element_1D(ib - double(Se.NB(1)) + 1)
                continue;
            end

            % -------------------------------------------------------------
            % PORT-NOTE [SCATTER_Bamp_gpu]
            %   call SCATTER_Bamp_gpu(isc)
            % fills isc%rhotw(1:ng) for transition (n,k) <- (m,k-q) with os(1)=ib
            % -------------------------------------------------------------
            rho = complex(randn(ng, 1, 'double'), randn(ng, 1, 'double')) * double(0.1);

            % isc%os in Fortran: (ob, ok, os, o_sp) — use toy f at (ib, ik, isp)
            ik  = 1;
            isp = 1;
            f_occ = double(E_f(ib, ik, isp));

            DP_Sc = cohsex_accumulate_toy(M_eff, rho, f_occ, spin_occ, DP_Sc);
        end
    end

    % Fortran: SE(ID)%SF(i_bas,1) += cmplx(DP_Sc(1)) + cmplx(DP_Sc(2))
    net = DP_Sc(1) + DP_Sc(2);
    Se.SF(i_bas, 1) = Se.SF(i_bas, 1) + complex(double(net), double(0));
    fprintf('  i_bas=%d: SEX=%.6g COH=%.6g net=%.6g\n', i_bas, DP_Sc(1), DP_Sc(2), net);
end

fprintf('\ndemo_cohsex_minimal: SF (real parts accumulated) = \n');
disp(real(Se.SF(:, 1)));

%% --- Local helpers (same algebra as GW_COHSEX_*.F) ------------------------
function Mout = cohsex_residuals_toy(M, gamp)
    % GW_COHSEX_residuals: Mout(ig1,ig2) = M(ig1,ig2) * gamp(ig1,ig2)
    Mout = M .* gamp;
end

function DP_Sc = cohsex_accumulate_toy(M, rho, f_occ, spin_occ, DP_Sc)
    % GW_COHSEX_accumulate (full row/col block here: rows 1:ng, cols 1:ng)
    rho_c = conj(rho);
    % local_rhotw = M * conj(rho)  (Fortran GEMV 'N')
    local_rhotw = M * rho_c;
    % V_dot_V on full slice -> use same contraction as row-local dot
    pre_factor = double(real(sum(conj(rho) .* local_rhotw)));

    pi_ = pi;
    DP_Sc(1) = DP_Sc(1) - (4 / spin_occ) * pi_ * pre_factor * f_occ;  % SEX
    DP_Sc(2) = DP_Sc(2) + 2 * pi_ * pre_factor;                        % COH
end

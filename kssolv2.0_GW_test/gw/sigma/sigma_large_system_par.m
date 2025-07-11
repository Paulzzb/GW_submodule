function sig = sigma_large_system_par(eps, sig, sys, options, syms)
%% Initialize essential values only
ryd = 13.6056923;
nbands = sig.nbnd;
ndiag_max = sig.ndiag_max;
nspin = sys.nspin;
sig.qpt = options.kpts;
sig.nkn = sys.nkpts;

%% Initialize grid and vectors
sigrid = Ggrid(sys, 4 * sys.ecutwfc);
gvec = Gvector(sigrid,sys);
gr = fullbz(sys, syms, true);
fact = 1/(gr.nf * sys.vol);
coulfact = 8 * pi * fact;

%% Precompute frequently used values
wfc_cutoff = sys.ecutwfc * 2;
no_symmetries_q_grid = sig.no_symmetries_q_grid;
coul_cutoff = sig.coul_cutoff;

%% Initialize output arrays directly
sig.cor = zeros(ndiag_max, sig.nkn, nspin);
sig.sig = zeros(ndiag_max, sig.nkn, nspin);

%% Precompute k-point dependent values
ekin = zeros(gvec.ng, sig.nkn);
sig.isrtx = zeros(gvec.ng, sig.nkn);
sig.nmtx = zeros(1, sig.nkn);
sig.mtx = cell(sig.nkn, 1);

for ik = 1:sig.nkn
    rk = sig.qpt(ik, :);
    [ekin(:,ik), sig.isrtx(:,ik)] = sortrx(rk, gvec.ng, gvec.mill, sys);
    sig.nmtx(ik) = gcutoff(gvec.ng, ekin(:,ik), sig.isrtx(:,ik), eps.cutoff);
    sig.mtx{ik} = gvec.mill(sig.isrtx(1:sig.nmtx(ik), ik), :);
end

%% Precompute FBZ values
for ik = 1 : gr.nf
    rk = gr.f(ik, :);
    [ekin(:,ik), fbz.isrtx(:,ik)] = sortrx(rk, gvec.ng, gvec.mill, sys);
    fbz.nmtx(:,ik) = gcutoff(gvec.ng, ekin(:,ik), fbz.isrtx(:,ik), wfc_cutoff);
    fbz.mtx{ik} = gvec.mill(fbz.isrtx(1:fbz.nmtx(ik), ik), :);
    fbz.nmtx_cutoff(:,ik) = gcutoff(gvec.ng, ekin(:,ik), fbz.isrtx(:,ik), eps.cutoff);
    fbz.mtx_cutoff{ik} = gvec.mill(fbz.isrtx(1:fbz.nmtx_cutoff(ik), ik), :);
    
    % gmap for eps_inv
    isorti = zeros(gvec.ng, 1);
    for i = 1 : gvec.ng
        isorti(sig.isrtx(i, gr.indr(ik))) = i;
    end
    fbz.isorti(:, ik) = isorti;
end

fprintf('Starting sigma calculation loop over spins and bands...\n');

%% Main calculation loop - parallelized over bands (in)
for ispin = 1:nspin
    fprintf('Processing spin %d of %d...\n', ispin, nspin);
    
    % Temporary storage for parallel results
    cor_temp = zeros(ndiag_max, sig.nkn);
    sig_temp = zeros(ndiag_max, sig.nkn);
    
    parfor in = 1:ndiag_max
        fprintf('Processing band %d of %d...\n', in, ndiag_max);
        
        % Local copies for parfor
        local_cor = zeros(1, sig.nkn);
        local_sig = zeros(1, sig.nkn);
        
        for ik = 1:sig.nkn
            rk = sig.qpt(ik, :);
            syms_rk = subgrp(rk, syms);
            [nrk, neq, indrk] = irrbz(syms_rk, gr);
            wfnk = genwf(rk, gr, gvec, syms, sys, options, wfc_cutoff);
            
            if no_symmetries_q_grid
                nrk = gr.nf;
                indrk = (1:nrk)';
                neq = ones(nrk, 1);
            end
            
            % Preallocate accumulators
            asx_sum = 0;
            ax_sum = 0;
            ach_sum = 0;
            achx_sum = 0;

            for iq = 1:nrk
                ik_tmp = gr.indr(indrk(iq));
                itran = gr.itran(indrk(iq));
                qk = gr.r(ik_tmp,:) * syms.mtrx{itran,:};
                [~, kgq] = krange(qk, 1e-9);
                indt = gmap(gvec, syms, sig.nmtx(:,ik_tmp), itran, kgq, fbz.isrtx(:,indrk(iq)) ,fbz.isorti(:, indrk(iq)), sys);
                eps_inv = eps.inv{ik_tmp}(indt, indt);
                n_cutoff = fbz.nmtx_cutoff(indrk(iq));
                
                if ~no_symmetries_q_grid
                    [nstar, ~, ~] = rqstar(syms_rk, gr.f(indrk(iq), :));
                    if nstar ~= neq(iq)
                        error('nstar of kpoint %d mismatch', gr.f(indrk(iq), :))
                    end
                end
                
                % Compute Coulomb potential
                coulg = getvcoul(fbz.nmtx(indrk(iq)), fbz.isrtx(:, indrk(iq)), ...
                                ekin(:, indrk(iq)), coul_cutoff, 1);
                coulg_nocut = fact * coulg;
                coulg_cutoff = coulg_nocut(1:n_cutoff);
                
                % Compute eps_inv_I and related matrices
                I = eye(n_cutoff);
                eps_inv_I = eps_inv - I;
                eps_inv_I_coul = eps_inv_I .* coulg_cutoff';
                
                % Get wavefunction of k-q
                rkq = rk - gr.f(indrk(iq), :);
                wfnkq = genwf(rkq, gr, gvec, syms, sys, options, wfc_cutoff);
                occ_kq = get_occ(options, wfnkq.ikq, ispin);
                
                % Initialize local accumulators
                asx_loc = 0;
                ax_loc = 0;
                ach_loc = 0;
                
                % Vectorized band sum where possible
                for nn = 1:nbands
                    aqs_nn = getm_sigma(in, nn, wfnkq, wfnk, indrk(iq), ispin, fbz, sys);
                    aqs_cutoff = aqs_nn(1:n_cutoff);
                    aqs_product = aqs_cutoff .* aqs_cutoff';
                    aqs_eps_coul = aqs_product .* eps_inv_I_coul;
                    
                    if occ_kq(nn) > 0
                        asx_loc = asx_loc - occ_kq(nn) * aqs_eps_coul;
                        ax_loc = ax_loc - occ_kq(nn) * abs(aqs_nn).^2 .* coulg_nocut;
                    end
                    ach_loc = ach_loc + aqs_eps_coul;
                end
                
                % Exact static CH correlation if needed
                if sig.exact_static_ch
                    if indrk(iq) == 1
                        aqsch_in = getm_sigma(in, in, wfnkq, wfnk, 1, ispin, fbz, sys);
                    end
                    G_q = fbz.mtx_cutoff{indrk(iq)};
                    ncouls = fbz.nmtx_cutoff(1, indrk(iq));
                    
                    [ig, igp] = meshgrid(1:ncouls, 1:ncouls);
                    linearidx = sub2ind([ncouls ncouls], ig, igp);
                    gpp = zeros(ncouls * ncouls, 3);
                    gpp(linearidx, :) = bsxfun(@minus, G_q(ig, :), G_q(igp, :));
                    igpp = findvector(gpp, gvec);
                    igpp = fbz.isorti(igpp, 1);
                    
                    valid_indices = (igpp >= 1 & igpp <= fbz.nmtx(1));
                    aqsch_tmp = zeros(ncouls);
                    aqsch_tmp(valid_indices) = aqsch_in(igpp(valid_indices));
                    achx_loc = aqsch_tmp .* eps_inv_I_coul;
                    achx_sum = achx_sum + sum(achx_loc,"all") * neq(iq);
                end
                
                % Accumulate results
                asx_sum = asx_sum + sum(asx_loc,"all") * neq(iq);
                ax_sum = ax_sum + sum(ax_loc,"all") * neq(iq);
                ach_sum = ach_sum + sum(ach_loc,"all") * neq(iq);
            end
            
            % Store local results
            if sig.exact_static_ch
                local_cor(ik) = real(asx_sum + 0.5*achx_sum) * ryd;
                local_sig(ik) = real(asx_sum + ax_sum + 0.5*achx_sum) * ryd;
            else
                local_cor(ik) = real(asx_sum + 0.5*ach_sum) * ryd;
                local_sig(ik) = real(asx_sum + ax_sum + 0.5*ach_sum) * ryd;
            end
        end
        
        % Store results for this band
        cor_temp(in, :) = local_cor;
        sig_temp(in, :) = local_sig;
    end
    
    % Combine results from parallel workers
    sig.cor(:,:,ispin) = cor_temp;
    sig.sig(:,:,ispin) = sig_temp;
end

%% Final calculations
emf = ryd * options.ev;
sig = quasi_energy(nspin, ndiag_max, emf, sys.vxc, sig);

fprintf('Calculation completed.\n');
end
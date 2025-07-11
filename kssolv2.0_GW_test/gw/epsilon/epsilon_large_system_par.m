function eps = epsilon_large_system_par(sys, options, syms, eps)
%% Initialize some value
nvbands = eps.nv;
ncbands = eps.nc;
nbands = eps.nbnd;
nspin = sys.nspin;
nspinor = sys.nspinor;
wfc_cutoff = sys.ecutwfc * 2; % Ha --> Ry
fprintf('System parameters: nvbands = %d, ncbands = %d, nbands = %d, nspin = %d, nspinor = %d\n', nvbands, ncbands, nbands, nspin, nspinor);

sigrid = Ggrid(sys, 4*sys.ecutwfc);
gvec = Gvector(sigrid, sys);

pol.qpt = options.kpts;

%%
gr = fullbz(sys, syms, true);

ekin = zeros(gvec.ng, sys.nkpts);
for iq = 1:sys.nkpts
    qq = pol.qpt(iq,:);
    [ekin(:,iq), pol.isrtx(:,iq)] = sortrx(qq, gvec.ng, gvec.mill, sys);
    pol.nmtx(:,iq) = gcutoff(gvec.ng, ekin(:,iq), pol.isrtx(:,iq), eps.cutoff);
    pol.mtx{:, iq} = gvec.mill(pol.isrtx(1:pol.nmtx(iq), iq), :);
    
    %% get fftsize for calculating M matrix
    eps_box_min = zeros([1 3]);
    eps_box_max = zeros([1 3]);
    [eps_box_min, eps_box_max] = get_gvecs_bounds(pol.mtx{:, iq}, eps_box_min, eps_box_max);
    pol.fftgrid{:, iq} = min((options.wfn_fftgrid + eps_box_max - eps_box_min), options.fftgrid);
end

eps_tmp = cell([sys.nkpts 1]);
eps_inv = cell([sys.nkpts 1]);
chi0 = cell([sys.nkpts sys.nspin]);

for iq = 1:sys.nkpts
    fprintf('Calculating chi0 and epsilon for k-point %d/%d...\n', iq, sys.nkpts);
    qq = pol.qpt(iq,:);
    syms_qq = subgrp(qq, syms);
    fact = 4 / (gr.nf * sys.vol * nspin * nspinor);
    [nrq, neq, indrk] = irrbz(syms_qq, gr);
    
    % Initialize chi0 for current q-point
    current_size = pol.nmtx(1, iq);
    chi0_sum = zeros(current_size, current_size, nspin);
    
    % Precompute symmetry mappings outside parfor
    isorti_all = cell(nrq, 1);
    indt_all = cell(nrq, 1);
    rk_all = cell(nrq, 1);
    nstar_all = zeros(nrq, 1);
    for ik = 1:nrq
        rk = gr.f(indrk(ik),:);
        rk_all{ik} = rk;
        [nstar, indst, ~] = rqstar(syms_qq, rk);
        nstar_all(ik) = nstar;
        
        isorti = zeros(gvec.ng, 1);
        for i = 1:gvec.ng
            isorti(pol.isrtx(i, iq), 1) = i;
        end
        isorti_all{ik} = isorti;
        
        indt = cell(nstar, 1);
        for it = 1:nstar
            itran = syms_qq.indsub(indst(it));
            kgq = -syms_qq.kgzero(indst(it),:);
            indt{it} = gmap(gvec, syms, pol.nmtx(:,iq), itran, kgq, pol.isrtx(:,iq), isorti, sys);
        end
        indt_all{ik} = indt;
    end
    
    parfor ik = 1:nrq
        fprintf('Processing irreducible k-point %d/%d in k-point %d/%d...\n', ik, nrq, iq, sys.nkpts);
        rk = rk_all{ik};
        nstar = nstar_all(ik);
        indt = indt_all{ik};
        
        %% get wavefunction of k+q and k
        rkq = rk + qq;
        wfnkq = genwf(rkq, gr, gvec, syms, sys, options, wfc_cutoff);
        wfnk = genwf(rk, gr, gvec, syms, sys, options, wfc_cutoff);
        
        temp_chi0 = zeros(current_size, current_size, nspin);
        
        for ispin = 1:nspin
            occ_vkq = get_occ(options, wfnkq.ikq, ispin);
            no_v = sum(occ_vkq > 0);
            occ_ck = get_occ(options, wfnk.ikq, ispin);
            no_c = sum(occ_ck > 0) + 1;
            
            for iv = 1:no_v
                for ic = no_c:nbands
                    % Compute M matrix element
                    m_element = getm_epsilon(iv, ic, wfnkq, wfnk, iq, ispin, pol, sys, options);
                    temp_chi0(:,:,ispin) = temp_chi0(:,:,ispin) - conj(m_element) * m_element.';
                    
                    % Handle degenerate k-points
                    if (nstar > 1)
                        for it = 2:nstar
                            m_element_degen = m_element(indt{it});
                            temp_chi0(:,:,ispin) = temp_chi0(:,:,ispin) - conj(m_element_degen) * m_element_degen.';
                        end
                    end
                end
            end
        end
        
        % Accumulate using reduction (MATLAB automatically handles this)
        chi0_sum = chi0_sum + temp_chi0;
    end
    
    % Finalize chi0 for this q-point
    for ispin = 1:nspin
        chi0{iq, ispin} = chi0_sum(:,:,ispin) * fact;
    end
    
    if (nspin == 2)
        chi0{iq,1} = chi0{iq, 1} + chi0{iq, 2};
    end
    
    % Compute epsilon
    coulg = getvcoul(pol.nmtx(1, iq), pol.isrtx(:, iq), ekin(:, iq), eps.coul_cutoff, 0);
    eps_tmp{iq} = eye(current_size) - coulg .* chi0{iq, 1};
    eps_inv{iq} = inv(eps_tmp{iq});
end

eps.inv = eps_inv;
eps.mtx = pol.mtx;
eps.nmtx = pol.nmtx;
fprintf('Calculation completed. Results stored in eps.inv and eps.mtx.\n');
end
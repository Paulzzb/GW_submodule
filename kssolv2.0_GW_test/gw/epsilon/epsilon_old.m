function eps = epsilon(sys, options, syms, eps)
%% Initialize some value
nvbands = eps.nv;
ncbands = eps.nc;
nbands = eps.nbnd;
nspin = sys.nspin;
nspinor = sys.nspinor;
sigrid = Ggrid(sys, 4*sys.ecutwfc);
gvec = Gvector(sigrid,sys);
pol.qpt = options.kpts;

%%
gr = fullbz(sys, syms, true);
gme = cell([nvbands ncbands gr.nf sys.nkpts sys.nspin]);

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


for iq = 1 : sys.nkpts
    qq = pol.qpt(iq,:);
    syms_qq = subgrp(qq, syms);
    fact= 4 / (gr.nf * sys.vol * nspin * nspinor);
    [nrq, neq, indrk] = irrbz(syms_qq,gr);
    
    for ik = 1 : nrq
        rk = gr.f(indrk(ik),:);
        [nstar, indst, rqs] = rqstar(syms_qq, rk);
        if (nstar ~= neq(ik))
            error('nstar of kpoint %d mismatch', rk)
        end
        for it = 1:nstar
            itran = syms_qq.indsub(indst(it));
            kgq = -syms_qq.kgzero(indst(it),:);
            for i = 1:gvec.ng
                isorti(pol.isrtx(i, iq), 1) = i;
            end
            indt{it, ik} = gmap(gvec, syms, pol.nmtx(:,iq), itran, kgq, pol.isrtx(:,iq) ,isorti, sys);
        end
        
        %% get wavefunction of k+q
        rkq = rk + qq;
        wfnkq = genwf(rkq, gr, gvec, pol, syms, sys, options, eps.cutoff);
        wfnk = genwf(rk, gr, gvec, pol, syms, sys, options, eps.cutoff);
        for ispin = 1:nspin
            occ_vkq = get_occ(options, wfnkq.ikq, ispin);
            no_v = sum(occ_vkq > 0);
            occ_ck = get_occ(options, wfnk.ikq, ispin);
            no_c = sum(occ_ck > 0) + 1;
            for iv = 1:no_v
                for ic = no_c:nbands
                    gme{iv,ic,indrk(ik),iq,ispin} = getm_epsilon(iv, ic, wfnkq, wfnk, iq, ispin, pol, sys, options);
                    %% If there are degenerate kpoints in k space, find it and assign it
                    if (nstar > 1)
                        for it = 2:nstar
                            k_degenerate = rqs(it, :);
                            [~, ik_degenerate] = ismember(k_degenerate, gr.f, 'rows');
                            gme{iv,ic,ik_degenerate,iq,ispin} = gme{iv,ic,indrk(ik),iq,ispin}(indt{it, ik});
                        end
                    end
                end
            end
        end
    end
end

eps_tmp = cell([sys.nkpts 1]);
eps_inv = cell([sys.nkpts 1]);
chi0 = cell([sys.nkpts sys.nspin]);
for iq = 1:sys.nkpts
    for ispin = 1:nspin
        gme_q = gme(:, :, :, iq, ispin);
        chi0_tmp = cellfun(@(x) conj(x) * x.' * fact, gme_q, 'UniformOutput', false); % Use the cellfun function to operate on each element in a cell array
        chi0{iq,ispin} = -sum(cat(3,chi0_tmp{:}),3); % The negative sign comes from the extra negative sign of taking the conjugate of the energy denominator.
    end
    if (nspin == 2)
        chi0{iq,1} = chi0{iq, 1} + chi0{iq, 2}; % Sum the spin components
    end
    coulg = getvcoul(pol.nmtx(1, iq), pol.isrtx(:, iq), ekin(:, iq), eps.coul_cutoff, 0);
    eps_tmp{iq} = eye(pol.nmtx(1, iq));
    eps_tmp{iq} = eps_tmp{iq}- coulg .* chi0{iq, 1};
    eps_inv{iq} = inv(eps_tmp{iq});
end

eps.inv = eps_inv;
eps.mtx = pol.mtx;
end
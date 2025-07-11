function sig = sigma(eps, sig, sys, options, syms)
%% Initialize some value
ryd = 13.6056923;
nbands = sig.nbnd;
ndiag_min = sig.ndiag_min;
ndiag_max = sig.ndiag_max;
wfc_cutoff = sys.ecutwfc * 2; % Ha --> Ry
nspin = sys.nspin;
nspinor = sys.nspinor;

ndiag = ndiag_max - ndiag_min;
aqs = cell(nbands, nspin);
aqsch = cell(nbands, nspin);
asx = zeros([ndiag sys.nkpts nspin]);
ax = zeros([ndiag sys.nkpts nspin]);
ach = zeros([ndiag sys.nkpts nspin]);
achx = zeros([ndiag sys.nkpts nspin]);
sigrid = Ggrid(sys, 4 * sys.ecutwfc);
gvec = Gvector(sigrid,sys);
coul_cutoff = sig.coul_cutoff;
no_symmetries_q_grid = sig.no_symmetries_q_grid;
sig.qpt = options.kpts;
sig.nkn = sys.nkpts;

%%
gr = fullbz(sys, syms, true);
fact = 1/(gr.nf * sys.vol);
coulfact = 8 * pi * fact;
eps_inv_fbz = cell([gr.nf 1]);

for ik = 1 : sig.nkn
    rk = sig.qpt(ik, :);
    [ekin(:,ik), sig.isrtx(:,ik)] = sortrx(rk, gvec.ng, gvec.mill, sys);
    sig.nmtx(:,ik) = gcutoff(gvec.ng, ekin(:,ik), sig.isrtx(:,ik), eps.cutoff);
    sig.mtx{:, ik} = gvec.mill(sig.isrtx(1:sig.nmtx(ik), ik), :); % TODO: the G vector grid should be the same with the one in eps
end

for ik = 1 : gr.nf
    rk = gr.f(ik, :);
    [ekin(:,ik), fbz.isrtx(:,ik)] = sortrx(rk, gvec.ng, gvec.mill, sys);
    fbz.nmtx(:,ik) = gcutoff(gvec.ng, ekin(:,ik), fbz.isrtx(:,ik), wfc_cutoff);
    fbz.mtx{:, ik} = gvec.mill(fbz.isrtx(1:fbz.nmtx(ik), ik), :);
    fbz.nmtx_cutoff(:,ik) = gcutoff(gvec.ng, ekin(:,ik), fbz.isrtx(:,ik), eps.cutoff);
    fbz.mtx_cutoff{:, ik} = gvec.mill(fbz.isrtx(1:fbz.nmtx_cutoff(ik), ik), :);
    %% gmap for eps_inv
    
    itran = gr.itran(ik);
    qk = gr.r(gr.indr(ik),:) * syms.mtrx{itran,:};
    [~, kgq] = krange(qk, 1e-9);
    for i = 1 : gvec.ng
        isorti(sig.isrtx(i, gr.indr(ik)), 1) = i;
    end
    fbz.isorti(:, ik) = isorti;
    indt = gmap(gvec, syms, sig.nmtx(:,gr.indr(ik)), itran, kgq, fbz.isrtx(:,ik) ,isorti, sys);
    eps_inv_fbz{ik} = eps.inv{gr.indr(ik)}(indt, indt);
end

fprintf('Starting sigma calculation loop over spins and bands...\n');

for ispin = 1 : nspin
    fprintf('Processing spin %d of %d...\n', ispin, nspin);
    for in = ndiag_min : ndiag_max
        fprintf('Processing band %d ...\n', in);
        for ik = 1 : sig.nkn
            fprintf('Processing k-point %d of %d...\n', ik, sig.nkn);
            asxtemp = 0;
            axtemp = 0;
            achtemp = 0;
            achxtemp = 0;
            rk = sig.qpt(ik, :);
            syms_rk = subgrp(rk, syms);
            [nrk, neq, indrk] = irrbz(syms_rk, gr);
            wfnk = genwf(rk, gr, gvec, syms, sys, options, wfc_cutoff);
            if no_symmetries_q_grid
                nrk = gr.nf;
                indrk = (1:nrk);
                neq = ones(1, nrk);
            end
            for iq = 1 : nrk
                asx_loc = 0;
                ax_loc = 0;
                ach_loc = 0;
                achx_loc = 0;
                n_cutoff = fbz.nmtx_cutoff(1, indrk(iq));
                qq = gr.f(indrk(iq), :);
                eps_inv = eps_inv_fbz{indrk(iq)};
                if ~no_symmetries_q_grid
                    [nstar, indst, rqs] = rqstar(syms_rk, qq);
                    if (nstar ~= neq(iq))
                        error('nstar of kpoint %d mismatch', qq)
                    end
                end
                
                %%
                I = eye(fbz.nmtx_cutoff(indrk(iq)));
                coulg = getvcoul(fbz.nmtx(1, indrk(iq)), fbz.isrtx(:, indrk(iq)), ekin(:, indrk(iq)), coul_cutoff, 1);
                coulg_nocut = fact * coulg;
                coulg_cutoff = coulg_nocut(1 : n_cutoff, 1);
                eps_inv_I = eps_inv - I;
                eps_inv_I_coul = eps_inv_I .* coulg_cutoff';
                %% get wavefunction of k-q
                rkq = rk - qq;
                wfnkq = genwf(rkq, gr, gvec, syms, sys, options, wfc_cutoff);
                %% Sum over band nn
                occ_kq = get_occ(options, wfnkq.ikq, ispin);
                for nn = 1 : nbands
                    aqs{nn, ispin} = getm_sigma(in, nn, wfnkq, wfnk, indrk(iq), ispin, fbz, sys);
                    aqs_nocut = aqs{nn, ispin};
                    aqs_cutoff = aqs{nn, ispin}(1 : n_cutoff, 1);
                    aqs_product = aqs_cutoff.* aqs_cutoff';
                    aqs_eps_coul = aqs_product .* eps_inv_I_coul;
                    %% Calculate SX and X
                    if occ_kq(nn) > 0
                        asx_loc = asx_loc - occ_kq(nn) * aqs_eps_coul; % this is actually asx - x;
                        ax_loc = ax_loc - occ_kq(nn) * abs(aqs_nocut).^2 .* coulg_nocut;
                    end
                    %% Calculate CH without exact ch correlation
                    ach_loc = ach_loc + aqs_eps_coul;
                end
                
                %% Calculate CH with exact ch correlation
                if sig.exact_static_ch
                    if (indrk(iq) == 1)
                        %aqsch{in,ispin} = getm_sigma(in, in, wfnkq, wfnk, 1, ispin, fbz, sys);
                        aqsch{in, ispin} = aqs{in, ispin};
                    end
                    G_q = fbz.mtx_cutoff{1, indrk(iq)};
                    ncouls = fbz.nmtx_cutoff(1, indrk(iq));
%                     index_tmp = zeros(ncouls);
                    aqsch_tmp = zeros(ncouls);
                    isorti = fbz.isorti(:, 1);
                    
                    [ig, igp] = meshgrid(1:ncouls, 1:ncouls);
                    linearidx = sub2ind([ncouls ncouls], ig, igp);
                    
                    gpp = zeros(ncouls * ncouls, 3);
                    gpp(linearidx, :) = bsxfun(@minus, G_q(ig, :), G_q(igp, :));
                    %gpp(linearidx, :) = G_q(ig, :) - G_q(igp, :);
                    igpp = findvector(gpp, gvec);
                    igpp = isorti(igpp);
                    
                    valid_indices = (igpp >= 1 & igpp <= fbz.nmtx(1, 1));
%                     index_tmp(valid_indices) = igpp(valid_indices);
                    aqsch_tmp(valid_indices) = aqsch{in, ispin}(igpp(valid_indices));
                    aqsch_tmp(~valid_indices) = 0;
                    
%                     TOL_Small = 1e-6;
%                     eps_inv_I(abs(eps_inv_I) < TOL_Small) = 0;
%                     eps_inv_I_coul = eps_inv_I .* coulg_cutoff'; %TODO: Is it needed to remove small part in eps_inv?
                    
                    achx_loc = aqsch_tmp .* eps_inv_I_coul;
                    achxtemp = achxtemp + sum(achx_loc,"all") * neq(iq);
                end
                asxtemp = asxtemp + sum(asx_loc,"all") * neq(iq);
                axtemp = axtemp + sum(ax_loc,"all") * neq(iq);
                achtemp = achtemp + sum(ach_loc,"all") * neq(iq);
            end
            n_index = in - ndiag_min + 1;
            asx(n_index,ik,ispin) = sum(asxtemp,"all");
            ax(n_index,ik,ispin) = sum(axtemp,"all");
            ach(n_index,ik,ispin) = 0.5 * sum(achtemp,"all");
            if sig.exact_static_ch
                achx(n_index,ik,ispin) = 0.5 * sum(achxtemp,"all");
            end
        end
    end
end

fprintf('Finalizing calculations...\n');

if sig.exact_static_ch
    sig.cor = real(asx + achx) * ryd;
    sig.sig = real(asx + ax + achx) * ryd;
else
    sig.cor = real(asx + ach) * ryd;
    sig.sig = real(asx + ax + ach) * ryd;
end

emf = ryd * options.ev;
sig = quasi_energy(nspin, ndiag_min, ndiag_max, emf, sys.vxc, sig);

fprintf('Calculation completed.\n');
end
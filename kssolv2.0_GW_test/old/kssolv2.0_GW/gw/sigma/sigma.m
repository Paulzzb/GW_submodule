function sig = sigma(eps, sig, sys, options, syms)
%% Initialize some value
ryd = 13.6056923;
nbands = sig.nbnd;
ndiag_max = sig.ndiag_max;
aqs = cell([nbands 1]);
aqsch = cell([nbands 1]);
asx = zeros([ndiag_max sys.nkpts]);
ax = zeros([ndiag_max sys.nkpts]);
ach = zeros([ndiag_max sys.nkpts]);
achx = zeros([ndiag_max sys.nkpts]);
nspin = sys.nspin;
nspinor = sys.nspinor;
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
    fbz.nmtx(:,ik) = gcutoff(gvec.ng, ekin(:,ik), fbz.isrtx(:,ik), eps.cutoff);
    fbz.mtx{:, ik} = gvec.mill(fbz.isrtx(1:fbz.nmtx(ik), ik), :);
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

for ispin = 1 : nspin
    for in = 1 : ndiag_max
        for ik = 1 : sig.nkn
            asxtemp = 0;
            axtemp = 0;
            achtemp = 0;
            achxtemp = 0;
            rk = sig.qpt(ik, :);
            syms_rk = subgrp(rk, syms);
            [nrk, neq, indrk] = irrbz(syms_rk, gr);
            wfnk = genwf(rk, gr, gvec, sig, syms, sys, options, eps.cutoff);
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
                qq = gr.f(indrk(iq), :);
                if ~no_symmetries_q_grid
                    [nstar, indst, rqs] = rqstar(syms_rk, qq);
                    if (nstar ~= neq(iq))
                        error('nstar of kpoint %d mismatch', qq)
                    end
                end
                
                %%
                I = eye(fbz.nmtx(indrk(iq)));
                coulg = getvcoul(fbz.nmtx(1, indrk(iq)), fbz.isrtx(:, indrk(iq)), ekin(:, indrk(iq)), coul_cutoff, 1);
                coulg = fact * coulg;
                %% get wavefunction of k-q
                rkq = rk - qq;
                wfnkq = genwf(rkq, gr, gvec, sig, syms, sys, options, eps.cutoff);
                %% Sum over band nn
                for nn = 1 : nbands
                    aqs{nn,ispin} = getm_sigma(in, nn, wfnkq, wfnk, indrk(iq), ispin, fbz, sys);
                    %% Calculate SX and X
                    occ_kq = get_occ(options, wfnkq.ikq, ispin);
                    asx_loc = asx_loc - occ_kq(nn) * aqs{nn,ispin} .* aqs{nn,ispin}' .* eps_inv_fbz{indrk(iq)} * coulg;
                    ax_loc = ax_loc - occ_kq(nn) * aqs{nn,ispin}.* aqs{nn,ispin}'.* I * coulg;
                    %ax_loc=ax_loc-occ_kq(nn)*abs(diag(aqs{nn,ispin})).^2.*coulg;
                    %% Calculate CH without exact ch correlation
                    ach_loc = ach_loc + 0.5 * aqs{nn,ispin}.* aqs{nn,ispin}'.* (eps_inv_fbz{indrk(iq)} -I) * coulg;
                end
                
                %% Calculate CH with exact ch correlation
                if sig.exact_static_ch
                    if (indrk(iq) == 1)
                        aqsch{in,ispin} = getm_sigma(in, in, wfnkq, wfnk, 1, ispin, sig, sys);
                    end
                    G_q = fbz.mtx{1, indrk(iq)};
                    ncouls = fbz.nmtx(1, indrk(iq));
                    index_tmp = eye(ncouls);
                    aqsch_tmp = eye(ncouls);
                    isorti = fbz.isorti(:, 1);
                    
                    for ig = 1:ncouls
                        for igp = 1:ncouls
                            linearidx(ig, igp) = sub2ind([ncouls ncouls], ig, igp);
                            gpp(linearidx(ig, igp), :) = G_q(ig,:) - G_q(igp,:);
                        end
                    end
                    [igpp, ~] = findvector2(gpp, gvec);
                    igpp = isorti(igpp);
                    for ig = 1:ncouls
                        for igp = 1:ncouls
                            igpp_tmp = igpp(linearidx(ig, igp), :);
                            if igpp_tmp >=1 && igpp_tmp <= sig.nmtx(1, 1)
                                index_tmp(ig, igp) = igpp_tmp;
                                aqsch_tmp(ig, igp) = aqsch{in, ispin}(igpp_tmp);
                            else
                                index_tmp(ig, igp) = 0;
                                aqsch_tmp(ig, igp) = 0;
                            end
                        end
                    end
                    achx_loc = 0.5 * aqsch_tmp .* (eps_inv_fbz{indrk(iq)} - I) * coulg;
                    achxtemp = achxtemp + sum(achx_loc,"all") * neq(iq);
                end
                asxtemp = asxtemp + sum(asx_loc,"all") * neq(iq);
                axtemp = axtemp + sum(ax_loc,"all") * neq(iq);
                achtemp = achtemp + sum(ach_loc,"all") * neq(iq);
            end
            asx(in,ik,ispin) = ryd * sum(asxtemp,"all");
            ax(in,ik,ispin) = ryd * sum(axtemp,"all");
            ach(in,ik,ispin) = ryd * sum(achtemp,"all");
            if sig.exact_static_ch
                achx(in,ik,ispin) = ryd * sum(achxtemp,"all");
            end
        end
    end
end
if sig.exact_static_ch
    sig.cor = real(asx - ax + achx);
    sig.sig = real(asx + achx);
else
    sig.cor = real(asx - ax + ach);
    sig.sig = real(asx + ach);
end
emf = ryd * options.ev;
sig = quasi_energy(nspin, ndiag_max, emf, sys.vxc, sig);
end
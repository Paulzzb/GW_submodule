function bsemat = kernel(eps, xct, sys, options, syms)
%% Initialize some values
nv = xct.nv;
nc = xct.nc;
xct.iscreen = 1;

xctrid = Ggrid(sys, 4 * sys.ecutwfc);
gvec = Gvector(xctrid, sys);
%%
gr = fullbz(sys, syms, true);
xct.qpt = gr.f;

for ik = 1 : sys.nkpts
    rk = options.kpts(ik, :);
    [ekin(:,ik), ibz.isrtx(:,ik)] = sortrx(rk, gvec.ng, gvec.mill, sys);
    ibz.nmtx(:,ik) = gcutoff(gvec.ng, ekin(:,ik), ibz.isrtx(:,ik), eps.cutoff);
    ibz.mtx{:, ik} = gvec.mill(ibz.isrtx(1:ibz.nmtx(ik), ik), :); % TODO: the G vector grid should be the same with the one in eps
end

for ik = 1 : gr.nf
    for ikp = 1 : gr.nf
        rk = xct.qpt(ik, :);
        rkp = xct.qpt(ikp, :);
        
        q.r = rk - rkp;
        q.length = norm(q.r);
        [ikq, itqq, kgqq] = find_kpt_match(gr, syms, q.r);
        [ekin, isrtc_kq] = sortrx(q.r, gvec.ng, gvec.mill, sys);
        q.nmtx = gcutoff(gvec.ng, ekin, isrtc_kq, eps.cutoff);
        q.mtx = gvec.mill(isrtc_kq(1 : q.nmtx), :);
        q.isort = isrtc_kq;
        for i = 1 : gvec.ng
            isorti(ibz.isrtx(i,ikq), 1) = i;
        end
        indt = gmap(gvec, syms, q.nmtx, itqq, kgqq, q.isort, isorti, sys);
        eps_inv_q = eps.inv{ikq}(indt, indt);
        coulg_q = getvcoul(q.nmtx, q.isort, ekin, xct.coul_cutoff, 0);
        coulg_q = coulg_q / (8 * pi);
        
        wfnk = genwf(rk, gr, gvec, ibz, syms, sys, options, eps.cutoff);
        wfnkp = genwf(rkp, gr, gvec, ibz, syms, sys, options, eps.cutoff);
        
        for ispin = 1 : sys.nspin
            for ispin_p = 1 : sys.nspin
                occ_vk_s = get_occ(options, wfnk.ikq, ispin);
                no_vk_s = sum(occ_vk_s > 0);
                occ_ck_s = get_occ(options, wfnk.ikq, ispin);
                no_ck_s = sum(occ_ck_s > 0) + 1;
                occ_vkp_s = get_occ(options, wfnkp.ikq, ispin);
                no_vkp_s = sum(occ_vkp_s > 0);
                occ_ckp_s = get_occ(options, wfnkp.ikq, ispin);
                no_ckp_s = sum(occ_ckp_s > 0) + 1;
                
                occ_vk_sp = get_occ(options, wfnk.ikq, ispin_p);
                no_vk_sp = sum(occ_vk_sp > 0);
                occ_ck_sp = get_occ(options, wfnk.ikq, ispin_p);
                no_ck_sp = sum(occ_ck_sp > 0) + 1;
                occ_vkp_sp = get_occ(options, wfnkp.ikq, ispin_p);
                no_vkp_sp = sum(occ_vkp_sp > 0);
                occ_ckp_sp = get_occ(options, wfnkp.ikq, ispin_p);
                no_ckp_sp = sum(occ_ckp_sp > 0) + 1;
                
                ind = 0;
                
                for iv = 1 : nv
                    actul_v_s = no_vk_s - nv + iv;
                    actul_v_sp = no_vk_sp - nv + iv;
                    for ic = 1 : nc
                        actul_c_s = no_ck_s + ic - 1;
                        actul_c_sp = no_ck_sp + ic - 1;
                        for ivp = 1 : nv
                            actul_vp_s = no_vkp_s - nv + ivp;
                            actul_vp_sp = no_vkp_sp - nv + ivp;
                            for icp = 1 : nc
                                actul_cp_s = no_ckp_s + icp - 1;
                                actul_cp_sp = no_ckp_sp + icp - 1;
                                
                                % Compute matrix elements: <ck|exp(i(k-kp-G0+G).r)|ckp> -> mccp
                                m_ck_ckp = getm_kernel(actul_c_s, actul_cp_s, wfnk, wfnkp, q, ispin, sys);
                                % Compute matrix element: <vk|exp(i(k-kp-G0+G).r)|vkp> -> mvvp
                                m_vkp_vk = getm_kernel(actul_vp_sp, actul_v_sp, wfnkp, wfnk, q, ispin_p, sys);
                                % Compute matrix elements: <vk|e^{i*G.r}|ck> -> mvc
                                m_vk_ck = getm_kernel(actul_v_s, actul_c_s, wfnk, wfnk, q, ispin, sys);
                                % Compute matrix elements: <vkp|e^{i*G.r}|ckp> -> mvpcp
                                m_vkp_ckp = getm_kernel(actul_vp_sp, actul_cp_sp, wfnkp, wfnkp, q, ispin_p, sys);
                                
                                w_mat = eps_inv_q .* coulg_q';
                                if (xct.iscreen == 1)
                                    w_mat(1:q.nmtx, 1) = w_mat(1:q.nmtx, 1) * q.length;
                                    w_mat(1, 1:q.nmtx) = w_mat(1, 1:q.nmtx) * q.length;
                                    w_mat(1, 1) = 1;
                                end
                                k_direct_head = m_ck_ckp(1) .* w_mat(1, 1) .* m_vkp_vk(1);
                                k_direct_wings = m_ck_ckp(1) .* w_mat(1, 2 : q.nmtx) .* m_vkp_vk(2 : q.nmtx) + ...,
                                    m_ck_ckp(2 : q.nmtx)' .* w_mat(2 : q.nmtx, 1) .* m_vkp_vk(1);
                                k_direct_body = m_ck_ckp(2 : q.nmtx)' .* w_mat(2 : q.nmtx, 2 : q.nmtx) .* m_vkp_vk(2 : q.nmtx);
                                k_direct = m_ck_ckp' .* w_mat .* m_vkp_vk;
                                k_exchange = conj(m_vk_ck(2 : q.nmtx)) .* coulg_q(2 : q.nmtx) .* m_vkp_ckp(2 : q.nmtx);
                                
                                ind = ind + 1;
                                % bsehead(iv, ic, ik, ivp, icp, ikp, ispin, ispin_p) = sum(k_direct_head, 'all');
                                % bsewing(iv, ic, ik, ivp, icp, ikp, ispin, ispin_p) = sum(k_direct_wings, 'all');
                                % bsebody(iv, ic, ik, ivp, icp, ikp, ispin, ispin_p) = sum(k_direct_body, 'all');
                                % bsedirect(iv, ic, ik, ivp, icp, ikp, ispin, ispin_p) = sum(k_direct, 'all');
                                % bsex(iv, ic, ik, ivp, icp, ikp, ispin, ispin_p) = sum(k_exchange, 'all');
                                bsehead(ind, ispin, ispin_p) = sum(k_direct_head, 'all');
                                bsewing(ind, ispin, ispin_p) = sum(k_direct_wings, 'all');
                                bsebody(ind, ispin, ispin_p) = sum(k_direct_body, 'all');
                                bsedirect(ind, ispin, ispin_p) = sum(k_direct, 'all');
                                bsex(ind, ispin, ispin_p) = sum(k_exchange, 'all');
                                
                            end
                        end
                    end
                end
            end
        end
    end
end
bsemat.head = bsehead;
bsemat.wing = bsewing;
bsemat.body = bsebody;
bsemat.x = bsex;
end
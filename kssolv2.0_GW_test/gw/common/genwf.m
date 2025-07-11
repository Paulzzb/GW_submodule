function [wfnkq] = genwf(rkq, gr, gvec, syms, sys, options, cutoff)

[ikrkq, itqq, kgqq] = find_kpt_match(gr, syms, rkq);
rkmatch = gr.r(ikrkq, :);
[ekin, isrtc_kq] = sortrx(rkq, gvec.ng, gvec.mill, sys);
nmtx = gcutoff(gvec.ng, ekin, isrtc_kq, cutoff);

kadd = g2fft_index(options.mill{1, ikrkq}, sys); % The index of G vector in WFN file
isortc = gvec.index_vec(kadd);
for i = 1:nmtx
    isorti(isortc(i, 1), 1) = i;
end
ind = gmap(gvec, syms, nmtx, itqq, kgqq, isrtc_kq ,isorti, sys); % The index connect G of rkq and the one in ikrkq WFN
%% assignment
wfnkq.ikq = ikrkq;
wfnkq.mill = gvec.mill(isrtc_kq(1:nmtx), :);
wfnkq.ind = ind;
wfnkq.psi = cell([sys.nspin sys.nspinor]);
for ispin = 1 : sys.nspin
    for ispinor = 1 : sys.nspinor
        wfn_ispinor = options.X0.wavefuncell{ikrkq, ispin}.psi( 1 + (ispinor - 1) * nmtx : ispinor * nmtx, :);
        wfnkq.psi{ispin, ispinor} = wfn_ispinor(ind, :);
    end
end

%% In spinor case, we must rotate spinors according to spinor rotation matrix umtrx
if sys.nspinor == 2
    if itqq ~= 1
        umtrx = susymmetries(sys.bvec, syms.mtrx{itqq}, itqq);
        for ispin = 1 : sys.nspin
            for iband = 1 : sys.nbnd
                wfn_iband = [wfnkq.psi{ispin, 1}(:, iband) wfnkq.psi{ispin, 2}(:, iband)] * umtrx;
                wfnkq.psi{ispin, 1}(:, iband) = wfn_iband(:, 1);
                wfnkq.psi{ispin, 2}(:, iband) = wfn_iband(:, 2);
            end
            wfn_ispin = [wfnkq.psi{ispin, 1} ;wfnkq.psi{ispin, 2}];
            checknorm(wfn_ispin, ikrkq, ispin); % we check if the norm differs appreciably from unity. there is no longer a need to further normalize the vector
        end
    end
end

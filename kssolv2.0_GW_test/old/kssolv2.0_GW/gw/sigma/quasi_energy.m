function sig = quasi_energy(nspin, ndiag_max, emf, vxc, sig)
%% Symmetrize matrix elements
for ik = 1:sig.nkn
    for ispin = 1:nspin
        nsubs = 1;
        ndeg(nsubs) = 1;
        for ii = 2:ndiag_max
            dek = emf(ii, ik, ispin) - emf(ii-1, ik, ispin);
            tolerance = 1e-5;  % 设置接近1的容忍度
            if abs(dek)<tolerance
                ndeg(nsubs)=ndeg(nsubs)+1;
            else
                nsubs = nsubs+1;
                ndeg(nsubs) = 1;
            end
        end
        istop = 0;
        for ii = 1:nsubs
            istart = istop+1;
            istop = istart+ndeg(ii)-1;
            acor = 0;
            asig = 0;
            for jj = istart:istop
                acor = acor + sig.cor(jj, ik, ispin);
                asig = asig + sig.sig(jj, ik, ispin);
            end
            fact = 1/ndeg(ii);
            for jj = istart:istop
                sig.cor(jj, ik, ispin) = acor*fact;
                sig.sig(jj, ik, ispin) = asig*fact;
            end
        end
    end
    %% Get quasiparticle energy
    [~, ind] = ismember(sig.qpt(ik, :), vxc.kpoints, 'rows');
    if (nspin == 2)
        sig.vxc(:, ik, 1) = vxc.value(1:2:2*ndiag_max, ind);
        sig.vxc(:, ik, 2) = vxc.value(2:2:2*ndiag_max, ind);
    else
        sig.vxc(:, ik) = vxc.value(1:ndiag_max, ind);
    end
end
sig.eqp0 = emf(1:ndiag_max, :, :) - sig.vxc + sig.sig;
end
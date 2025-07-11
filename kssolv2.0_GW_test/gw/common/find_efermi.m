function efermi = find_efermi(rfermi, efermi_input, nband, minband, search_efermi, update_efermi, sys, options)
% rfermi: whether or not adjust according to efermi_input
% efermi_input: relative or absolute Fermi level
ev = options.ev;
ifmax = options.ifmax;
if ~search_efermi && update_efermi
    error('BUG: cannot call find_efermi with should_update but not should_search')
end
if nband < 1 || nband > size(ev, 1)
    error('find_efermi: nband out of bounds')
end
if all(ifmax(1 : sys.nkpts, 1 : sys.nspin) <= 0)
    error('All k-points have no occupied bands')
end

if search_efermi
    %% get middle energy, not differentiate between spins
    for ik = 1 : sys.nkpts
        for ispin = 1 : sys.nspin
            vbm(ik, ispin) = max(ev(minband : ifmax(ik, ispin), ik, ispin));
            cbm(ik, ispin) = min(ev(ifmax(ik, ispin) + 1 : nband, ik, ispin));
        end
    end
    emiddle = (max(max(vbm)) + min(min(cbm))) / 2;
    %% check. Question: any inconsistency caused by the spin energy level?
    for ik = 1 : sys.nkpts
        for ispin = 1 : sys.nspin
            if ifmax(ik, ispin) > 0
                if any(ev(minband : ifmax(ik, ispin), ik, ispin) > emiddle)
                    error('there is a valence state above the middle energy')
                end
            elseif ifmax(ik, ispin) + 1 <= nband
                if any(ev(ifmax(ik, ispin) + 1 : ifmax(ik, ispin), ik, ispin) > emiddle)
                    error('there is a conduction state below the middle energy')
                end
            end
        end
    end
    %% adjust
    if rfermi
        efermi_tmp = emiddle + efermi_input;
    else
        efermi_tmp = emiddle;
    end
    if update_efermi
        efermi = efermi_tmp;
    end
end
%TODO: reset occupations(reset ifmax if Fermi level was moved by input
%file, or using Fermi level from other wfns). Is it really necessary?

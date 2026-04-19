function wf = get_wf_at_bz(ib, ikbz, ispin, k_data)
ikibz = int32(k_data.bz2ibz(ikbz));
irot = int32(k_data.bz2rot(ikbz));
isc = int32([ib, ikibz, irot, ispin]);
wf = double(wave_functions.WF_apply_symm(isc));
end

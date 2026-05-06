function out = get_coeff(id, param)


isdf_data = isdf.get(id);
symm_data = symmetry.get;
r_lat_data = lattice.manager('r_lat', 'get');
fft_data = FFT.manager('get');

bs = isdf_data.bundle_struct;

ib = param.isc(1);
ikibz = param.isc(2);
ikrot = param.isc(3);
ispin = param.isc(4);
ob = param.iscp(1);
ikpibz = param.iscp(2);
ikprot = param.iscp(3);

iGo = param.iGo;
iqibz = param.iqibz;
iqrot = param.iqrot;


if ~strcmp(isdf_data.interp_scheme, 'coarse')
  iscoarse = false;
else
  iscoarse = true;
end


if ~iscoarse
  wf_ibz_in_bundle_a = bs.WF_bundle(:, ib, ikibz, ispin);
end
% wf_ibz = wf_data.c(:, ib, ikibz, ispin);
if ikrot > nsym / (is_t_rev + 1)
  ikrot_wf = ikrot -  nsym / (is_t_rev + 1);
  ikrot_wf = symm_data.inv_rot_index(ikrot_wf);
  conjflag1 = true;
else
  ikrot_wf = ikrot;
  ikrot_wf = symm_data.inv_rot_index(ikrot_wf);
  conjflag1 = false;
end
if ikprot > nsym / (is_t_rev + 1)
  ikprot_wf = ikprot -  nsym / (is_t_rev + 1);
  ikprot_wf = symm_data.inv_rot_index(ikprot_wf);
  conjflag2 = true;
else
  ikprot_wf = ikprot;
  ikprot_wf = symm_data.inv_rot_index(ikprot_wf);
  conjflag2 = false;
end

if ~strcmp(isdf_data.interp_scheme, 'coarse')
  indSq = bs.R_rot_in_bundle(:, iqrot);
  indices = bs.R_rot_in_bundle(:, ikrot_wf);
  wf_bz_in_bundle = wf_ibz_in_bundle_a(indices);
  if conjflag1
    wf_bz_in_bundle = conj(wf_bz_in_bundle);
  end
  wf_Sq_bundle = wf_bz_in_bundle(indSq);
  wf1_Sq_xalpha = wf_Sq_bundle(bs.sampling2bundle);
  
  %

  wf_ibz_in_bundle_a  = bs.WF_bundle(:, ob, ikpibz, ispin);
  indices = bs.R_rot_in_bundle(:, ikprot_wf);
  wf_bz_in_bundle = wf_ibz_in_bundle_a(indices);
  if conjflag2
    wf_bz_in_bundle = conj(wf_bz_in_bundle);
  end
  wf_Sq_bundle = wf_bz_in_bundle(indSq);
  wf2_Sq_xalpha = wf_Sq_bundle(bs.sampling2bundle);
else
  indSq = isdf_data.R_rot_coarse(:, iqrot);
  wf1_ikibz = isdf_data.coeff_seper(:, ib, ikibz, ispin);
  ind = isdf_data.R_rot_coarse(:, ikrot_wf);
  wf_1_ikbz = wf1_ikibz(ind);
  if conjflag1
    wf_1_ikbz = conj(wf_1_ikbz);
  end
  wf1_Sq_xalpha = wf_1_ikbz(indSq);

  wf2_ikpibz = isdf_data.coeff_seper(:, ob, ikpibz, ispin);
  ind = isdf_data.R_rot_coarse(:, ikprot_wf);
  wf_2_ikpibz = wf2_ikpibz(ind);
  if conjflag2
    wf_2_ikpibz = conj(wf_2_ikpibz);
  end
  wf2_Sq_xalpha = wf_2_ikpibz(indSq);
end

out = conj(wf1_Sq_xalpha) .* wf2_Sq_xalpha;

Go = single(r_lat_data.Ggrid_RLU(iGo, :));
inviqrot = symm_data.inv_rot_index(iqrot);
invSqGo = single(Go * symm_data.rot_mtrx_RLU_G(:, :, inviqrot));
phase_shift_coarse = exp(-2 * pi * 1i * (isdf_data.R_sampling_RLU ./ double(fft_data.fftgrid)) * invSqGo');
% phase_shift_coarse = exp(-2 * pi * 1i * SqR_scal * Go');
out = out .* phase_shift_coarse;

end
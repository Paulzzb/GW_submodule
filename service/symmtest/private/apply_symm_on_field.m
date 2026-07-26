function out = apply_symm_on_field(field_in, isym, symm_data, fft_data)
if isym == 1
  out = field_in;
  return;
end

inv_isym = int32(symm_data.inv_rot_index(isym));
ind = fft_data.R_rot(:, inv_isym);
out = field_in(ind);

nsym = double(symm_data.nsym);
is_t_rev = double(symm_data.is_t_rev);
if isym > nsym / (1 + is_t_rev)
  out = conj(out);
end
end

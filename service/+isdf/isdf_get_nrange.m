function [nrange1, nrange2] = isdf_get_nrange(id)
  % Get band index ranges (nrange1, nrange2) for ISDF type "nn" / "vn" / "vc".
  wf_data = wave_functions.get();
  isdf_data = isdf.get(id);
  type = isdf_data.desc;
  nb = wf_data.nb;

  if strcmp(type, "nn")
    nrange1 = 1:nb;
    nrange2 = 1:nb;
    return
  end

  system_data = system.get();
  nkibz = lattice.manager('k', 'get').nibz;
  nspin = system_data.nspin;
  nocc_max = 0;
  for ispin = 1:nspin
    for ikibz = 1:nkibz
      f_ib = system_data.f(:, ikibz, ispin);
      idx_last = find(f_ib(:) > 1e-5, 1, 'last');
      if ~isempty(idx_last)
        nocc_max = max(nocc_max, idx_last);
      end
    end
  end

  if strcmp(type, "vn")
    nrange1 = 1:nocc_max;
    nrange2 = 1:nb;
  elseif strcmp(type, "vc")
    nrange1 = 1:nocc_max;
    nrange2 = (nocc_max + 1):nb;
  else
    error('isdf_get_nrange:Type', 'Unknown isdf desc ''%s''.', type);
  end
end

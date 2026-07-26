function Vxc = load_vxc_from_dat(vxc_path, nb, nk, nspin)
%LOAD_VXC_FROM_DAT  Read QE vxc.dat into [nb x nk x nspin] (eV, as in the file).
%
%   Vxc = load_vxc_from_dat(vxc_path, nb, nk)       % nspin = 1
%   Vxc = load_vxc_from_dat(vxc_path, nb, nk, nspin)
%
%   File layout matches load_qe_from_folder.m (header row: kx ky kz nband ...).

  if nargin < 4 || isempty(nspin)
    nspin = 1;
  end
  if ~isfile(vxc_path)
    error('load_vxc_from_dat:NotFound', 'vxc.dat not found: %s', vxc_path);
  end

  fid = dlmread(vxc_path);
  [vxcrow, ~] = size(fid);
  if vxcrow < 2
    error('load_vxc_from_dat:Format', 'vxc.dat is too short: %s', vxc_path);
  end

  step = fid(1, 4);
  if step < 1 || ~isfinite(step)
    error('load_vxc_from_dat:Format', 'Invalid band count in vxc.dat header: %s', vxc_path);
  end

  vxc_kpts = [];
  vxc_vals = [];
  for i = 1:(step + 1):(vxcrow - step)
    vxc_kpts = [vxc_kpts; fid(i, 1:3)]; %#ok<AGROW>
    a = [];
    for j = i + 1:i + step
      a = [a, fid(j, 3)]; %#ok<AGROW>
    end
    vxc_vals = [vxc_vals, a]; %#ok<AGROW>
  end
  vxc_vals = vxc_vals.';

  nk_file = size(vxc_vals, 2);
  if nk_file ~= nk
    warning('load_vxc_from_dat:nk', ...
      'vxc.dat has %d k-points; requested nk=%d. Using min(nk).', nk_file, nk);
    nk_use = min(nk, nk_file);
  else
    nk_use = nk;
  end

  nbnd_vxc = floor(size(vxc_vals, 1) / nspin);
  nb_use = min(nb, nbnd_vxc);
  if nb_use < nb
    warning('load_vxc_from_dat:nb', ...
      'vxc.dat provides %d bands per spin; requested nb=%d. Padding with zeros.', ...
      nbnd_vxc, nb);
  end

  Vxc = zeros(nb, nk, nspin);
  for ik = 1:nk_use
    if nspin == 1
      Vxc(1:nb_use, ik, 1) = vxc_vals(1:nb_use, ik);
    else
      for ispin = 1:nspin
        vxc_spin = vxc_vals(ispin:nspin:end, ik);
        Vxc(1:nb_use, ik, ispin) = vxc_spin(1:nb_use);
      end
    end
  end
end

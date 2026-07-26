% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/26

function driver(data, config)
  if nargin < 1
    error('system.driver requires input data.');
  end
  if nargin < 2
    config = [];
  end

  if ~isfield(data, 'ev')
    error('system.driver requires data.ev for ground-state energies.');
  end

  if ~isfield(data, 'occupation')
    error('system.driver requires data.occupation for occupations.');
  end

  Eo = double(data.ev);
  f = double(data.occupation);
  if isfield(data, 'Vxc') && ~isempty(data.Vxc)
    Vxc = double(data.Vxc);
  elseif isfield(data, 'vxc') && ~isempty(data.vxc)
    Vxc = double(data.vxc);
  else
    Vxc = [];
  end
  if isempty(Vxc)
    Vxc = local_vxc_from_groundstate(config, size(Eo));
  end
  if isempty(Vxc)
    warning('No Vxc in groundstate data and no vxc.dat; system_data.Vxc set to zeros.');
    Vxc = zeros(size(Eo));
  end
  if ndims(Eo) == 2
    Eo = reshape(Eo, size(Eo, 1), size(Eo, 2), 1);
  end
  if ndims(f) == 2
    f = reshape(f, size(f, 1), size(f, 2), 1);
  end
  if isvector(Vxc)
    Vxc = reshape(Vxc, [], 1, 1);
  elseif ndims(Vxc) == 2
    Vxc = reshape(Vxc, size(Vxc, 1), size(Vxc, 2), 1);
  end

  if ~isequal(size(Eo), size(f))
    error('system.driver requires data.ev and data.occupation to have matching sizes.');
  end
  if ~isequal(size(Eo), size(Vxc))
    error('system.driver requires data.Vxc to have the same size as data.ev.');
  end

  [nb, nk, nspin] = size(Eo);
  system_data = system.base.system_m(nb, nk, nspin);
  system_data.Eo = Eo;
  system_data.Vxc = Vxc;
  system_data.f = f;
  system_data.qptype = "HF";
  system_data.assigned = true;
  system_data.allocated = true;
  system.save2mod(system_data);


  system.degeneracy_detect();

end

function Vxc = local_vxc_from_groundstate(config, eo_size)
  Vxc = [];
  if ~isstruct(config) || ~isfield(config, 'CONTROL') ...
      || ~isfield(config.CONTROL, 'groundstate_dir')
    return
  end
  vxc_path = fullfile(config.CONTROL.groundstate_dir, 'vxc.dat');
  if ~isfile(vxc_path)
    return
  end
  nb = eo_size(1);
  nk = eo_size(2);
  nspin = 1;
  if numel(eo_size) >= 3
    nspin = eo_size(3);
  end
  Vxc = load_vxc_from_dat(vxc_path, nb, nk, nspin);
end


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

  if ~isfield(data, 'ev')
    error('system.driver requires data.ev for ground-state energies.');
  end

  if ~isfield(data, 'occupation')
    error('system.driver requires data.occupation for occupations.');
  end

  Eo = double(data.ev);
  f = double(data.occupation);
  if ndims(Eo) == 2
    Eo = reshape(Eo, size(Eo, 1), size(Eo, 2), 1);
  end
  if ndims(f) == 2
    f = reshape(f, size(f, 1), size(f, 2), 1);
  end

  if ~isequal(size(Eo), size(f))
    error('system.driver requires data.ev and data.occupation to have matching sizes.');
  end

  [nb, nk, nspin] = size(Eo);
  system_data = system.base.system_m(nb, nk, nspin);
  system_data.Eo = Eo;
  system_data.f = f;
  system_data.qptype = "HF";
  system_data.assigned = true;
  system_data.allocated = true;
  system.save2mod(system_data);


  system.degeneracy_detect();
  



end


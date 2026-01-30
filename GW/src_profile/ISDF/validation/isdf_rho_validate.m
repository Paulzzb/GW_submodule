% 
% License-Identifier: GPL
% 
% Copyright (C) 2026
% 
% Authors (see AUTHORS file for details): ZZ 
% 
% Last modified: 2026/01/30 ZZ
function out = isdf_rho_validate(type, dbroot, GWinfo)
  %
  msg = "Validating ISDF for density calculations";
  QPlog(msg);
  %
  default_Constant = constant_map();
  nameConstants = fieldnames(default_Constant);
  for i = 1:numel(nameConstants)
    eval(sprintf('%s = %.16f;', nameConstants{i}, default_Constant.(nameConstants{i})));
  end
  %
  switch type
    case "vc"
      IID = 1; SID = "1";
    case "vs"
      IID = 2; SID = "2";
    case "ss"
      IID = 3; SID = "3";
  end
  %
  meta = db_read_meta(dbroot); 
  ind_xga = db_read(dbroot, IID, meta, "ind_xga");
  pga = db_read(dbroot, IID, meta, "pga");
  %
  psir = GWinfo.psir;
  vol = GWinfo.vol;
  occupation = GWinfo.occupation;
  nv = find(GWinfo.occupation > 1 - TOL_SMALL, 1, 'last');
  nr = size(psir, 1);
  tmp = "desc_type"+SID;
  Nisdf = meta.desc.(tmp).get("Nisdf");
  gvec = GWinfo.gvec;
  nfftgrid = gvec.nfftgridpts;
  n_start_end = meta.desc.(tmp).get("nlist");
  nlist = n_start_end(1):n_start_end(2);
  m_start_end = meta.desc.(tmp).get("mlist");
  mlist = m_start_end(1):m_start_end(2);
  %
  rhor = zeros(nr, 1);
  rhorisdf = zeros(nr, 1);
  c_rho = zeros(Nisdf, 1);
  %
  pgar = zeros(nr, Nisdf);
  for j = 1:Nisdf
    fftbox = put_into_fftbox(pga(:, j), gvec.idxnz, gvec.fftgrid);
    fftbox = do_FFT(fftbox, gvec.fftgrid, 1) * (nfftgrid / vol);
    pgar(:, j) = fftbox(:);
  end
  %
  for i = 1:nv
    rhor = rhor + occupation(i) * abs(psir(:, i)).^2;
    c_rho = abs(psir(ind_xga, i)).^2;
    rhorisdf = rhorisdf + occupation(i) * pgar * c_rho;
  end
  %
  norm_rho = norm(rhor, 2);
  norm_diff = norm(rhorisdf - rhor, 2);
  out = [norm_diff, norm_diff / norm_rho];
end

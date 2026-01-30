%
% License-Identifier: GPL
% 
% Copyright (C) 2026
% 
% Authors (see AUTHORS file for details): ZZ 
% 
% Last modified: 2026/01/30 ZZ
function out = isdf_hf_validate(type, dbroot, GWinfo)
  %
  msg = "Validating ISDF for HF calculations";
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
  hVh = db_read(dbroot, IID, meta, "hVh");
  %
  psir = GWinfo.psir;
  occupation = GWinfo.occupation;
  nv = find(GWinfo.occupation > 1 - TOL_SMALL, 1, 'last');
  ng = GWinfo.gvec.ng;
  Dcoul = spdiags(GWinfo.coulG, 0, ng, ng) * ry2ev;
  Dcoul(1,1) = GWinfo.coulG0 * ry2ev;
  vol = GWinfo.vol;
  %
  tmp = "desc_type"+SID;
  Nisdf = meta.desc.(tmp).get("Nisdf");
  n_start_end = meta.desc.(tmp).get("nlist");
  nlist = n_start_end(1):n_start_end(2);
  m_start_end = meta.desc.(tmp).get("mlist");
  mlist = m_start_end(1):m_start_end(2);
  btc = [min(m_start_end(1), n_start_end(1)),max(m_start_end(2), n_start_end(2))];
  %
  HF_dir  = zeros(btc(2) - btc(1) + 1, btc(2) - btc(1) + 1);
  HF_ISDF = zeros(btc(2) - btc(1) + 1, 1);
  %
  for ioper = 1:nv
    Mgvn = mtxel_sigma(ioper, GWinfo, btc(1):btc(2));
    Mgvn = conj(Mgvn);
    W1Mgvn = Dcoul * Mgvn;
    HF_dir = HF_dir + Mgvn' * W1Mgvn / vol;
  end
  HF_dir = diag(HF_dir);
  for i = 1:nv
    for j = btc(1):btc(2)
      c_rho = conj(psir(ind_xga, i)) .* psir(ind_xga, j);
      HF_ISDF(j) = HF_ISDF(j) + occupation(i) * c_rho' * hVh * c_rho;
    end
  end
  %
  diff = HF_dir - HF_ISDF;
  out = [max(abs(diff)), max(abs(diff./HF_dir))];
end

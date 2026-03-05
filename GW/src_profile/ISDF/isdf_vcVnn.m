% 
% License-Identifier: GPL
% 
% Copyright (C) 2026
% 
% Authors (see AUTHORS file for details): ZZ 
% 
% Last modified: 2026/02/02 ZZ
% function hVh = isdf_vcVnn(Phi, Psi, vcind, nnind, Dcoul, gvec, vol)
function hVh = isdf_vcVnn(dbroot, GWinfo, config)
% Calculate <p_mu | V | p_nu>, where mu <-- vc, and nu <-- nn
  cleanup = QPlog_push('ISDF-vcVnn');
  msg = sprintf('Log to be done');
  QPlog(msg, 2);
  %
  default_Constant = constant_map();
  nameConstants = fieldnames(default_Constant);
  for i = 1:numel(nameConstants)
    eval(sprintf('%s = %.16f;', nameConstants{i}, default_Constant.(nameConstants{i})));
  end
  %
  meta = db_read_meta(dbroot); 
  vcind = db_read(dbroot, 1, meta, "ind_xga");
  nnind = db_read(dbroot, 3, meta, "ind_xga");
  gvec = GWinfo.gvec;
  vol = GWinfo.vol;
  psir = GWinfo.psir;  ng = gvec.ng;
  Dcoul = spdiags(GWinfo.coulG, 0, ng, ng) * ry2ev;
  Dcoul(1,1) = GWinfo.coulG0 * ry2ev;
  %
  tmp = "desc_type"+"1";
  vcnrange = meta.desc.(tmp).get("nlist");
  vcnrange = vcnrange(1):vcnrange(2);
  vcmrange = meta.desc.(tmp).get("mlist");
  vcmrange = vcmrange(1):vcmrange(2);
  tmp = "desc_type"+"3";
  nnnrange = meta.desc.(tmp).get("nlist");
  nnnrange = nnnrange(1):nnnrange(2);
  nnmrange = meta.desc.(tmp).get("mlist");
  nnmrange = nnmrange(1):nnmrange(2);
  %
  Nblock = 64;
  % Nblock = 16*gcp('nocreate').NumWorkers;
  fftgrid = gvec.fftgrid;
  idxnz = gvec.idxnz;
  
  
  rkvc=length(vcind);
  phivc=psir(vcind, vcnrange);
  psivc=psir(vcind,vcmrange);
  C2=(phivc*phivc').*(psivc*psivc');
  if condest(C2) < 1e+12
    C2invvc = inv(C2);
  else
    C2invvc = pseudoinv(C2);
  end
  clear C2;
  
  rknn=length(nnind);
  phinn=psir(nnind,nnnrange);
  psinn=psir(nnind,nnmrange);
  C2=(phinn*phinn').*(psinn*psinn');
  if condest(C2) < 1e+12
    C2invnn = inv(C2);
  else
    C2invnn = pseudoinv(C2);
  end
  clear C2;
  %
  % -------------------------------------------------------------------
  % 1. Column slicing to calculate TMP = \F*(M C') 
  total_iter = ceil(rkvc / Nblock);
  est_total_time = -1;
  
  C1gvc = zeros(gvec.ng, rkvc);
  for i = 1:Nblock:rkvc
    iter_idx = ceil(i / Nblock);
  
    if iter_idx == 2
      startfirstiter = tic;
    end
  
    if i+Nblock < rkvc
      irange = i:i+Nblock-1;
    else
      irange = i:rkvc;
    end
  
    tmpr = (psir(:, vcnrange) * phivc(irange, :)') .* (psir(:, vcmrange) * psivc(irange, :)');
    for j = 0:length(irange)-1
    % parfor j = 0:length(irange)-1
      fftbox1 = reshape(tmpr(:, j+1), fftgrid);
      fftbox1 = do_FFT(fftbox1, fftgrid, 1) * vol;
      C1gvc(:, i+j) = get_from_fftbox(idxnz, fftbox1, fftgrid);
    end
  
    if iter_idx == 6 
      time_per_Nblock = toc(startfirstiter) / 4;
      est_total_time = time_per_Nblock * total_iter;
      msg = sprintf('Estimated total time: %.1f seconds', est_total_time);
      QPlog(msg, 2);
    end
  
    % Progress bar
    if est_total_time > 0
      elapsed_time = (iter_idx - 1) * time_per_Nblock;
      msg = sprintf('Progress: %3d%% | Elapsed: %.1fs / Estimated: %.1fs', ...
          round(100 * (iter_idx-1) / total_iter), elapsed_time, est_total_time);
      QPlog(msg, 2);
    else
      msg = sprintf('Progress: %3d%%', round(100 * (iter_idx-1) / total_iter));
      QPlog(msg, 2);
    end
  end
  %
  % -------------------------------------------------------------------
  % 1. Column slicing to calculate TMP = \F*(M C') 
  total_iter = ceil(rknn / Nblock);
  est_total_time = -1;
  
  C1gnn = zeros(gvec.ng, rknn);
  for i = 1:Nblock:rknn
    iter_idx = ceil(i / Nblock);
  
    if iter_idx == 2
      startfirstiter = tic;
    end
  
    if i+Nblock < rknn
      irange = i:i+Nblock-1;
    else
      irange = i:rknn;
    end
  
    tmpr = (psir(:, nnnrange) * phinn(irange, :)') .* (psir(:, nnmrange) * psinn(irange, :)');
    for j = 0:length(irange)-1
    % parfor j = 0:length(irange)-1
      fftbox1 = reshape(tmpr(:, j+1), fftgrid);
      fftbox1 = do_FFT(fftbox1, fftgrid, 1) * vol;
      C1gnn(:, i+j) = get_from_fftbox(idxnz, fftbox1, fftgrid);
    end
  
    if iter_idx == 6 
      time_per_Nblock = toc(startfirstiter) / 4;
      est_total_time = time_per_Nblock * total_iter;
      msg = sprintf('Estimated total time: %.1f seconds', est_total_time);
      QPlog(msg, 2);
    end
  
    % Progress bar
    if est_total_time > 0
      elapsed_time = (iter_idx - 1) * time_per_Nblock;
      msg = sprintf('Progress: %3d%% | Elapsed: %.1fs / Estimated: %.1fs', ...
          round(100 * (iter_idx-1) / total_iter), elapsed_time, est_total_time);
      QPlog(msg, 2);
    else
      msg = sprintf('Progress: %3d%%', round(100 * (iter_idx-1) / total_iter));
      QPlog(msg, 2);
    end
  end
  %
  % -------------------------------------------------------------------
  % 2. Calculate TMP2 = TMP'*Dcoul*TMP;
  C1vcVC1nn = zeros(rkvc, rknn);
  for i = 1:Nblock:rkvc
    iter_idx = ceil(i / Nblock);
    if iter_idx == 2
      startfirstiteri = tic;
    end
    if i+Nblock < rkvc
      irange = i:i+Nblock-1;
    else
      irange = i:rkvc;
    end
    %
    for j = 1:Nblock:rknn
      iter_jdx = ceil(j / Nblock);
      if iter_jdx == 2
        startfirstiterj = tic;
      end
  
      if j+Nblock < rknn
        jrange = j:j+Nblock-1;
      else
        jrange = j:rknn;
      end
      %
      C1vcVC1nn(irange, jrange) = C1gvc(:, irange)' * Dcoul * C1gnn(:, jrange) / vol;
    end
  end
  
  % -------------------------------------------------------------------
  % 4. Calculate hVh = TMP4' * TMP2 * TMP4
  hVh = zeros(rkvc, rknn);
  hVh = C2invvc * C1vcVC1nn * C2invnn;

end % function

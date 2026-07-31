% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/20 ZZ

function nm_Womega_nm = gw_fullfreq_cd_core_Gamma(config, n_start_end, m_start_end, omega_list, pattern)
%GW_FULLFREQ_CD_CORE_GAMMA  Single-(k,q) W matrix-element builder on service stack.
%
%   nm_Womega_nm = gw_fullfreq_cd_core_Gamma(config, [n1 n2], [m1 m2], omega_list, pattern)
%
% Computes <nm|W(q=Gamma;omega)|nm> for n in [n1,n2], m in [m1,m2], only
% for entries enabled by pattern (Nn x Nm x Nw).

default_Constant = constant_map();
nameConstants = fieldnames(default_Constant);
for i = 1:numel(nameConstants)
  eval(sprintf('%s = %.16f;', nameConstants{i}, default_Constant.(nameConstants{i})));
end

system_data = system.get();
k_data = lattice.manager('k', 'get');
q_data = lattice.manager('q', 'get');
r_lat_data = lattice.manager('r_lat', 'get');
coul_data = coulomb.get();

if k_data.nibz ~= 1 || q_data.nibz ~= 1 || q_data.nbz ~= 1
  error('gw_fullfreq_cd_core_Gamma:singleKQOnly', ...
    'Only double-k/double-q is supported (nkibz=%d, nqibz=%d, nqbz=%d).', ...
    k_data.nibz, q_data.nibz, q_data.nbz);
end

nstart = n_start_end(1);
nend = n_start_end(2);
mstart = m_start_end(1);
mend = m_start_end(2);
Nn = nend - nstart + 1;
Nm = mend - mstart + 1;
Nw = numel(omega_list);
mlist = mstart:mend;

if ~isequal(size(pattern), [Nn, Nm, Nw])
  error('gw_fullfreq_cd_core_Gamma:patternSize', ...
    'Pattern size must be [%d, %d, %d].', Nn, Nm, Nw);
end

ev = double(system_data.Eo(:, 1, 1)) * ry2ev;
focc = double(system_data.f(:, 1, 1));
nv = find(focc > 1 - TOL_SMALL, 1, 'last');
if isempty(nv)
  error('gw_fullfreq_cd_core_Gamma:occupation', 'Cannot determine nv from occupations.');
end

nb_total = numel(ev);
nsum = min(config.SYSTEM.number_bands_in_summation, nb_total);
if nsum <= nv
  error('gw_fullfreq_cd_core_Gamma:nsum', ...
    'Invalid summation range: nv=%d, nsum=%d (need nsum>nv).', nv, nsum);
end

nv_list = 1:nv;
nc_list = (nv + 1):nsum;
nc = numel(nc_list);

iqibz = q_data.bz2ibz(1);
iqrot = q_data.bz2rot(1);
iGo = r_lat_data.qindx_S(1, 1, 2);

vcoul_q = double(coul_data.vcoul(:, iqibz));
if iqibz == 1
  vcoul_q(1) = double(coul_data.vcoul0);
end
ng = numel(vcoul_q);
Dcoul = spdiags(vcoul_q(:), 0, ng, ng);
Dcoul = Dcoul * ry2ev;

% ISDF branch (preliminary): use vc/nn slots from service +isdf pool.
if isfield(config, 'ISDF') && isfield(config.ISDF, 'isisdf') && config.ISDF.isisdf
  use_isdf = true;
  [id_vc, ~, id_nn] = isdf.cohsex_resolve_ids(config);
  if isempty(id_vc) || isempty(id_nn)
    warning('gw_fullfreq_cd_core_Gamma:isdfFallback', ...
      ['config.ISDF.isisdf=true but vc/nn slots are not both available. ', ...
       'Fallback to dense Gamma path.']);
    use_isdf = false;
  end

  if use_isdf
    isdf.set_nrange(id_vc, config.SYSTEM);
    isdf.set_nrange(id_nn, config.SYSTEM);

    % Input key is lowercased by read_input_param -> verify_w_isdf
    verify_W_isdf = false;
    if isfield(config, 'FULLFREQ') && isfield(config.FULLFREQ, 'verify_w_isdf')
      verify_W_isdf = logical(config.FULLFREQ.verify_w_isdf);
    end
    vc_data = isdf.get(id_vc);

    nm_Womega_nm = zeros(Nn, Nm, Nw);
    for ifreq = 1:Nw
      omega = omega_list(ifreq);
      % service/+isdf/gen_Kq currently uses energies in Ry.
      omega_ry = omega / ry2ev;
      Kq_ISDF = isdf.gen_Kq(id_vc, iqibz, omega_ry);
      flagherm = abs(real(omega)) < 1e-5;
      if ifreq == 1
        flagherm = false;
      end
      if flagherm
        Kq_ISDF = (Kq_ISDF + Kq_ISDF') / 2;
      end
      tildeWq_nn = isdf.gen_tildeWq(id_vc, iqibz, Kq_ISDF, id_nn, flagherm);

      % Optional: dense G-space W vs ISDF reconstruction (same idea as gw_cohsex_multi_k).
      % Uses Hartree/Ry Coulomb (no ry2ev on v) to match gen_Kq / gen_tildeWq. Only at ω≈0
      % (Hermitian K) the static χ matches gen_Kq_Gamma; set &FULLFREQ verify_w_isdf = .true.
      if verify_W_isdf && flagherm
        nkbz = k_data.nbz;
        ng_vc = size(vc_data.helperqG, 1);
        if ng_vc ~= ng
          warning('gw_fullfreq_cd_core_Gamma:verifyWng', ...
            'helperqG row count (%d) ~= ng from vcoul (%d); skip W verify.', ng_vc, ng);
        else
          nrangev = double(vc_data.nrange1);
          nrangec = double(vc_data.nrange2);
          ev_ry = double(system_data.Eo);
          spin_id = 1;
          scal = 4.0;
          chiq_G = zeros(ng, ng);
          nrangev_row = reshape(nrangev, 1, []);
          nrangec_row = reshape(nrangec, 1, []);
          ncb = numel(nrangec_row);
          for ikbz = 1:nkbz
            ikibz_k = k_data.bz2ibz(ikbz);
            ikrot_k = k_data.bz2rot(ikbz);
            ikq_bz = r_lat_data.qindx_X(iqibz, ikbz, 1);
            iGo_x = r_lat_data.qindx_X(iqibz, ikbz, 2);
            ikq_ibz = k_data.bz2ibz(ikq_bz);
            ikq_rot = k_data.bz2rot(ikq_bz);
            f_c = double(system_data.f(nrangec_row, ikq_ibz, spin_id));
            e_c = ev_ry(nrangec_row, ikq_ibz, spin_id);
            for iv = nrangev_row
              Mgvc_blk = zeros(ng, ncb);
              for jc_id = 1:ncb
                jc = nrangec_row(jc_id);
                pchk = struct();
                pchk.is = [iv, ikibz_k, ikrot_k, spin_id];
                pchk.os = [jc, ikq_ibz, ikq_rot, spin_id];
                pchk.qs = [iGo_x, iqibz, iqrot];
                Mgvc_blk(:, jc_id) = double(SCATTER_Bamp(pchk));
              end
              f_v = double(system_data.f(iv, ikibz_k, spin_id));
              e_v = ev_ry(iv, ikibz_k, spin_id);
              occ = f_v - f_c;
              
                % den = e_v - e_c;
              den = omega_ry - e_c - e_v;
              valid = (abs(occ) >= 1e-8) & (abs(den) >= 1e-12);
              if ~any(valid)
                continue;
              end
              coeff = occ(valid) ./ den(valid);
              Mgvc_valid = Mgvc_blk(:, valid);
              Mgvc_weighted = Mgvc_valid .* reshape(coeff, 1, []);
              chiq_G = chiq_G + scal * (Mgvc_weighted * Mgvc_valid');
            end
          end
          inveps = eye(ng) - ( diag(vcoul_q) * chiq_G);
          W_dense_solve = full(inveps \ diag(vcoul_q));
          W_v = W_dense_solve - diag(vcoul_q);
          helperqG_vc = double(vc_data.helperqG(:, :, iqibz));
          K_use = double(Kq_ISDF);
          W_isdf = -diag(vcoul_q) * helperqG_vc * (K_use \ eye(size(K_use))) * helperqG_vc' * diag(vcoul_q);
          denom_solve = max(norm(W_v, 'fro'), eps);
          diff_solve = norm(W_v - W_isdf, 'fro');
          fprintf(['[gw_fullfreq_cd_core_Gamma verifyW] ifreq=%d iqibz=%d ng=%d Nisdf_nn=%d\n' ...
            '  ||W_dense-W_isdf||_F=%.6e (rel=%.6e)\n'], ...
            ifreq, iqibz, ng, size(tildeWq_nn, 1), diff_solve, diff_solve / denom_solve);
        end
      end

      for n = nstart:nend
        in = n - nstart + 1;
        row_pattern = pattern(in, :, ifreq);
        indm = find(row_pattern > 0);
        if isempty(indm)
          continue;
        end
        mlisttmp = mlist(indm);

        rho_mat = zeros(size(tildeWq_nn, 1), numel(mlisttmp));
        param = struct();
        param.is = [n, 1, 1, 1];
        param.qs = [iGo, iqibz, iqrot];
        for im = 1:numel(mlisttmp)
          m = mlisttmp(im);
          param.os = [m, 1, 1, 1];
          rho_mat(:, im) = isdf.get_rho_xalpha(id_nn, param);
        end

        Wrho = tildeWq_nn * rho_mat;
        out_list = sum(conj(rho_mat) .* Wrho, 1).';
        nm_Womega_nm(in, indm, ifreq) = out_list;
      end
    end
    nm_Womega_nm = nm_Womega_nm * ry2ev;
    return;
  end
end

Mvc_cache = cell(nv, 1);
Eden_cache = cell(nv, 1);

for iv = nv_list
  Mgvc = zeros(ng, nc);
  p = struct();
  p.is = [iv, 1, 1, 1];
  p.qs = [iGo, iqibz, iqrot];
  for ic = 1:nc
    cb = nc_list(ic);
    p.os = [cb, 1, 1, 1];
    Mgvc(:, ic) = conj(double(SCATTER_Bamp(p)));
  end
  Mvc_cache{iv} = Mgvc;
  Eden_cache{iv} = ev(iv) - ev(nc_list);
end

eta = 0.025;
nm_Womega_nm = zeros(Nn, Nm, Nw);

for ifreq = 1:Nw
  omega = omega_list(ifreq);
  if abs(real(omega)) < 1e-8
    ishermW = true;
  else
    ishermW = false;
  end

  chi_acc = zeros(ng, ng);
  for iv = nv_list
    Mgvc = Mvc_cache{iv};
    Eden = Eden_cache{iv};
    edenDR = (-1.0 ./ (omega - Eden - 1i * eta) + 1.0 ./ (omega + Eden + 1i * eta));
    chi_acc = chi_acc + 2.0 * Mgvc * (edenDR .* Mgvc') ;
  end

  if ishermW
    epsKernel = (Dcoul ) \ eye(ng) - chi_acc ;
    epsKernel = tril(epsKernel, -1) + tril(epsKernel, -1)' + diag(real(diag(epsKernel)));
    [L, D] = ldl(epsKernel);
    rsqrtD = diag(sqrt(diag(D)).^(-1));
  else
    A = eye(ng) - Dcoul * chi_acc;
    % Avoid explicit inverse: W = (I - A^{-1})Dcoul / vol.
    W = (Dcoul - (A \ Dcoul));
  end

  for n = nstart:nend
    in = n - nstart + 1;
    row_pattern = pattern(in, :, ifreq);
    indm = find(row_pattern > 0);
    if isempty(indm)
      continue;
    end
    mlisttmp = mlist(indm);

    Mnm = zeros(ng, numel(mlisttmp));
    p = struct();
    p.is = [n, 1, 1, 1];
    p.qs = [iGo, iqibz, iqrot];
    for im = 1:numel(mlisttmp)
      m = mlisttmp(im);
      p.os = [m, 1, 1, 1];
      Mnm(:, im) = conj(double(SCATTER_Bamp(p)));
    end

    if ishermW
      out_list = sum(Dcoul * abs(Mnm).^2, 1).';
      Mnm_r = rsqrtD * (L \ Mnm);
      out_list = out_list - sum(abs(Mnm_r).^2, 1).';
      nm_Womega_nm(in, indm, ifreq) = out_list;
    else
      WMnm = W * Mnm;
      out_list = sum(conj(Mnm) .* WMnm, 1).';
      nm_Womega_nm(in, indm, ifreq) = out_list;
    end
  end
end

end

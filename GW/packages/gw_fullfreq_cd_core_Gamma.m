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
d_lat_data = lattice.manager('d_lat', 'get');
coul_data = coulomb.get();

if k_data.nibz ~= 1 || q_data.nibz ~= 1 || q_data.nbz ~= 1
  error('gw_fullfreq_cd_core_Gamma:singleKQOnly', ...
    'Only single-k/single-q is supported (nkibz=%d, nqibz=%d, nqbz=%d).', ...
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

iqibz = q_data.bz2ibz(1);
iqrot = q_data.bz2rot(1);
iGo = r_lat_data.qindx_S(1, 1, 2);
vol = double(d_lat_data.DL_vol);

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

    nm_Womega_nm = zeros(Nn, Nm, Nw);
    for ifreq = 1:Nw
      omega = omega_list(ifreq);
      % service/+isdf/gen_Kq currently uses energies in Ry.
      omega_ry = omega / ry2ev;
      Kq_ISDF = isdf.gen_Kq(id_vc, iqibz, omega_ry);
      flagherm = abs(real(omega)) < 1e-5;
      tildeWq_nn = isdf.gen_tildeWq(id_vc, iqibz, Kq_ISDF, id_nn, flagherm);

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
    return;
  end
end

nv_list = 1:nv;
nc_list = (nv + 1):nsum;
nc = numel(nc_list);
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

eta = 0.0;
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
    chi_acc = chi_acc + 2.0 * Mgvc * (edenDR .* Mgvc') / vol;
  end

  if ishermW
    epsKernel = (Dcoul / vol) \ eye(ng) - chi_acc * vol;
    epsKernel = tril(epsKernel, -1) + tril(epsKernel, -1)' + diag(real(diag(epsKernel)));
    [L, D] = ldl(epsKernel);
    rsqrtD = diag(sqrt(diag(D)).^(-1));
  else
    A = eye(ng) - Dcoul * chi_acc;
    % Avoid explicit inverse: W = (I - A^{-1})Dcoul / vol.
    W = (Dcoul - (A \ Dcoul)) / vol;
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
      out_list = sum((Dcoul / vol) * abs(Mnm).^2, 1).';
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

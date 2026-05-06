function E = qp_cohsex(config)
%QP_COHSEX  Quasiparticle COHSEX route: Eqp = Ex + Esx_x + Ecoh (package path).

  QPlog('[COHSEX] Eqp = Ex + Esx_x + Ecoh (gw_x_k_packages + gw_cohsex_multi_k)...', 0);
  t0 = tic;

  Ex = gw_x_k_packages(config);
  [Esx_x, Ecoh] = gw_cohsex_multi_k(config);
  nik = size(Ex, 2);
  Eqp = zeros(size(Ex));
  for ik = 1:nik
    Eqp(:, ik) = Ex(:, ik) + Esx_x(:, ik) + Ecoh(:, ik);
    if ~isempty(Eqp(:, ik))
      msg = sprintf( ...
        '[COHSEX] ik=%d/%d  band1: Eqp=%.6f  Ex=%.6f  Esx_x=%.6f  Ecoh=%.6f', ...
        ik, nik, Eqp(1, ik), Ex(1, ik), Esx_x(1, ik), Ecoh(1, ik));
      QPlog(msg, 2);
    end
  end

  QPlog(sprintf('[COHSEX] Finished in %.2f s (nb=%d, nk=%d).', toc(t0), size(Eqp, 1), nik), 0);

  E.Eqp = Eqp;
  E.Ex = Ex;
  E.Esx_x = Esx_x;
  E.Ecoh = Ecoh;
end

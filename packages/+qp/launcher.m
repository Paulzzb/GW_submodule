% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/04 ZZ

function E = launcher(config)
%LAUNCHER  QP package dispatcher: Eqp = Ex + Esx_x + Ecoh by frequency_dependence.
%
%   E = qp.launcher(config)
%
% Routes:
%   frequency_dependence == -2  → gw.x_Gamma + gw.cohsex_Gamma
%   frequency_dependence == -1  → gw.x + gw.cohsex_multi_k
%   frequency_dependence ==  2  → full-frequency CD (res + int)

  msg = sprintf('[QP] Eqp = Ex + Esx_x + Ecoh (qp.launcher)...\n');
  output.msg('v0s', '%s', msg);
  t0 = tic;

  if config.FREQUENCY.frequency_dependence ~= -2
    Ex = gw.x(config);
  else
    Ex = gw.x_Gamma(config);
  end

  switch config.FREQUENCY.frequency_dependence
    case -2
      [Esx_x, Ecoh] = gw.cohsex_Gamma(config);
    case -1
      [Esx_x, Ecoh] = gw.cohsex_multi_k(config);
    case 2
      if ~isfield(config, 'freqinfo')
        config = generate_frequency(config);
      end
      Esx_x = gw.fullfreq_cd_res_Gamma(config);
      Ecoh = gw.fullfreq_cd_int_Gamma(config);
      ry2ev = constant_map().ry2ev;
      Esx_x = Esx_x / ry2ev;
      Ecoh = Ecoh / ry2ev;
    otherwise
      output.err('Invalid frequency dependence.');
  end
  nik = size(Ex, 2);
  Eqp = zeros(size(Ex));
  for ik = 1:nik
    Eqp(:, ik) = Ex(:, ik) + Esx_x(:, ik) + Ecoh(:, ik);
    if ~isempty(Eqp(:, ik))
      msg = sprintf( ...
        '[QP] ik=%d/%d  band1: Eqp=%.6f  Ex=%.6f  Esx_x=%.6f  Ecoh=%.6f\n', ...
        ik, nik, Eqp(1, ik), Ex(1, ik), Esx_x(1, ik), Ecoh(1, ik));
      output.msg('v2s', '%s', msg);
    end
  end

  msg = sprintf('[QP] Finished in %.2f s (nb=%d, nk=%d).\n', ...
    toc(t0), size(Eqp, 1), nik);
  output.msg('v0s', '%s', msg);

  E.Eqp = Eqp;
  E.Ex = Ex;
  E.Esx_x = Esx_x;
  E.Ecoh = Ecoh;
  E.Eqp0 = [];
  E.fout = '';

  try
    E.fout = qp.fout(E, config);
    msg = sprintf('[QP] Saved table to: %s\n', E.fout);
    output.msg('v0s', '%s', msg);
  catch ME
    msg = sprintf('[QP] Failed to write QP output file: %s\n', ME.message);
    output.msg('v0s', '%s', msg);
  end
end

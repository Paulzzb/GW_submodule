% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/06 ZZ

function E = launcher(config)
%LAUNCHER  QP package dispatcher: Eqp = Ex + Esx_x + Ecoh by frequency_dependence.
%
%   E = qp.launcher(config)
%
% Routes:
%   frequency_dependence == -2  → gw.x_Gamma + gw.cohsex_Gamma
%   frequency_dependence == -1  → gw.x + gw.cohsex_multi_k
%   frequency_dependence ==  2  → full-frequency CD (res + int)

  t0 = tic;

  if config.FREQUENCY.frequency_dependence == 2
    if ~isfield(config, 'freqinfo')
      config = generate_frequency(config);
    end
  end

  qp.summary(0, config);

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
      Esx_x = gw.fullfreq_cd_res_Gamma(config);
      Ecoh = gw.fullfreq_cd_int_Gamma(config);
    otherwise
      output.err('Invalid frequency dependence.');
  end

  nik = size(Ex, 2);
  Eqp = zeros(size(Ex));
  for ik = 1:nik
    Eqp(:, ik) = Ex(:, ik) + Esx_x(:, ik) + Ecoh(:, ik);
    if ~isempty(Eqp(:, ik))
      output.msg('v2s', ...
        'ik=%d/%d  band1: Eqp=%.6f  Ex=%.6f  Esx_x=%.6f  Ecoh=%.6f', ...
        ik, nik, Eqp(1, ik), Ex(1, ik), Esx_x(1, ik), Ecoh(1, ik));
    end
  end

  E.Eqp = Eqp;
  E.Ex = Ex;
  E.Esx_x = Esx_x;
  E.Ecoh = Ecoh;
  E.Eqp0 = [];
  E.fout = '';

  try
    E.fout = qp.fout(E, config);
  catch ME
    output.warn('Failed to write QP output file: %s', ME.message);
  end

  qp.summary(1, config, E, toc(t0));
end

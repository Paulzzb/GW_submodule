% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/04 ZZ

function fout = fout(E, config)
%FOUT  Write QP band table to a plain text file.
%
%   fout = qp.fout(E, config)
%
% Static (frequency_dependence == -1): one row per band (Emf, Eo, X, SX-X, CH, ...).
% Full frequency (frequency_dependence == 2): Re SX-X / Re CH / Re Sig / Re Eqp0 on
%   the first line; imaginary parts of SX-X, CH, Sig, and Eqp0 on the second line
%   (X stays on the first line only).

  def = constant_map();
  ry2ev = def.ry2ev;
  system_data = system.get();

  nbmin = config.SYSTEM.energy_band_index_min;
  nbmax = config.SYSTEM.energy_band_index_max;
  nband = nbmax - nbmin + 1;
  nk = size(E.Ex, 2);

  freq_dep = config.FREQUENCY.frequency_dependence;
  if freq_dep == 1
    freq_dep = 0;
  end

  Eo = double(system_data.Eo(nbmin:nbmax, 1:nk, 1)) * ry2ev;
  Vxc = double(system_data.Vxc(nbmin:nbmax, 1:nk, 1));
  Ex = double(E.Ex(:, 1:nk)) * ry2ev;
  Esx_x = double(E.Esx_x(:, 1:nk)) * ry2ev;
  Ecoh = double(E.Ecoh(:, 1:nk)) * ry2ev;
  Sig = Ex + Esx_x + Ecoh;
  Eqp0 = Eo + Sig - Vxc;
  Emf = Eo;

  out_dir = config.CONTROL.storage_dir;
  if isfield(config, 'CONTROL') && isfield(config.CONTROL, 'output_dir') ...
      && ~isempty(config.CONTROL.output_dir)
    out_dir = config.CONTROL.output_dir;
  end
  if ~exist(out_dir, 'dir')
    mkdir(out_dir);
  end

  fout = fullfile(out_dir, 'qp.dat');
  fid = fopen(fout, 'w');
  if fid == -1
    output.err('Cannot open output file: %s', fout);
  end

  cleaner = onCleanup(@() fclose(fid));

  switch freq_dep
    case 2
      fprintf(fid, '   n         Emf          Eo           X      Re SX-X       Re CH      Re Sig        Vxc     Re Eqp0\n');
      fprintf(fid, '                                     Im SX-X       Im CH      Im Sig                Im Eqp0\n');
      for ik = 1:nk
        for ib = 1:nband
          band_idx = nbmin + ib - 1;
          fprintf(fid, '%4d%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f \n', ...
            band_idx, Emf(ib, ik), Eo(ib, ik), real(Ex(ib, ik)), real(Esx_x(ib, ik)), ...
            real(Ecoh(ib, ik)), real(Sig(ib, ik)), Vxc(ib, ik), real(Eqp0(ib, ik)));
          fprintf(fid, '%40s%12.6f%12.6f%12.6f            %12.6f \n', ...
            '', imag(Esx_x(ib, ik)), imag(Ecoh(ib, ik)), imag(Sig(ib, ik)), imag(Eqp0(ib, ik)));
        end
      end
    otherwise
      fprintf(fid, '   n         Emf          Eo           X        SX-X          CH         Sig         Vxc        Eqp0\n');
      for ik = 1:nk
        for ib = 1:nband
          band_idx = nbmin + ib - 1;
          fprintf(fid, '%4d%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f\n', ...
            band_idx, Emf(ib, ik), Eo(ib, ik), Ex(ib, ik), Esx_x(ib, ik), ...
            Ecoh(ib, ik), Sig(ib, ik), Vxc(ib, ik), Eqp0(ib, ik));
        end
      end
  end
end

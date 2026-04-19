% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08

function [Esum2, EsumISDF2, DiffEsum2] = isdf_coarse_validate_energies(tildeVq, wf_on_coarse, R_rot_coarse, varargin)
% ISDF_COARSE_VALIDATE_ENERGIES  Compare direct Coulomb exchange-style energy vs ISDF tildeVq contraction.
%
% Mirrors the accumulation loops in gen_indices_coarse_test (SCATTER_Bamp vs c_rho' * tildeVq * c_rho).
% Optional trailing args:
%   One matrix R_sampling_RLU (Nwf x 3): same phase as gen_tildeVq, row-aligned with wf_on_coarse.
%   Or {R_coarse_RLU, fftgrid_i, fftgrid_c}: integer coarse box scaled by fftgrid_i./fftgrid_c (isdftest path).
%
% Progress and elapsed / estimated total time use timing.LIVE_timing (CPU time via timing.timing_string).

  wf_data = wave_functions.get();
  k_data = lattice.manager('k', 'get');
  r_lat_data = lattice.manager('r_lat', 'get');
  coul_data = coulomb.get();
  system_data = system.get();

  nibz = int32(k_data.nibz);
  nbz = int32(k_data.nbz);
  nb = int32(wf_data.nb);
  nspin = int32(wf_data.nspin);
  symm_data = symmetry.get();

  % Phase: (A) exactly one trailing arg = R_sampling_RLU (Nwf x 3), same as gen_tildeVq / driver path.
  %        (B) three or more trailing args = R_coarse_RLU, fftgrid_i, fftgrid_c (integer coarse box + scaling).
  use_coarse_G0_phase = (numel(varargin) >= 3);
  use_R_sampling = (numel(varargin) == 1 && ~isempty(varargin{1}));
  if use_R_sampling
    R_sampling_RLU = varargin{1};
  elseif use_coarse_G0_phase
    R_coarse_RLU = varargin{1};
    fftgrid_i = varargin{2};
    fftgrid_c = varargin{3};
  end

  Esum2 = 0.0;
  EsumISDF2 = 0.0;
  DiffEsum2 = 0.0;
  Esum = 0.0;
  DiffEsum = 0.0;

  if use_R_sampling && size(R_sampling_RLU, 1) ~= size(wf_on_coarse, 1)
    error('isdftest:isdf_coarse_validate_energies:BadR', ...
      'R_sampling_RLU rows (%d) must match wf_on_coarse (%d).', ...
      size(R_sampling_RLU, 1), size(wf_on_coarse, 1));
  end

  total_triples = max(1, double(nibz) * double(nspin) * double(nb));

  isdftest_ensure_timing_initialized();
  tm_live = timing.get();
  tm_live.live.nhash = int32(20);
  tm_live.live.live_report_min_seconds = 0;
  timing.save2mod(tm_live);

  timing.LIVE_timing('isdf validate (ik,ispin,ib)', total_triples);
  cleanup_live = onCleanup(@() timing.LIVE_timing());

  for ik = 1:nibz
    ikibz = ik;
    ikrot = 1;
    for ispin = 1:nspin
      for ib = 1:nb
        wf_1_c = wf_on_coarse(:, ib, ikibz, ispin);
        if ikrot ~= 1
          wf_1_c = isdf.isdf_apply_symm_on_coarse(wf_1_c, ikrot, symm_data, R_rot_coarse);
        end
        wf_1_c = conj(wf_1_c);
        for iqbz = 1:nbz
          iqibz = k_data.bz2ibz(iqbz, 1);
          iqrot = k_data.bz2rot(iqbz, 1);

          ikpbz = r_lat_data.qindx_S(ik, iqbz, 1);
          iGo = r_lat_data.qindx_S(ik, iqbz, 2);
          ikpibz = k_data.bz2ibz(ikpbz, 1);
          ikprot = k_data.bz2rot(ikpbz, 1);
          isc = [ib, ik, 1, ispin];


          vcoul_q = coul_data.vcoul(:, iqibz);
          if iqibz == 1
            vcoul_q(1) = coul_data.vcoul0;
          end


          for ob = 1:nb
            occupation = system_data.f(ob, ikpibz, ispin);
            if occupation < 1e-6
              continue;
            end
            wf_2_c = wf_on_coarse(:, ob, ikpibz, ispin);
            if ikprot ~= 1
              wf_2_c = isdf.isdf_apply_symm_on_coarse(wf_2_c, ikprot, symm_data, R_rot_coarse);
            end
            c_rho = wf_1_c .* wf_2_c;

            if use_R_sampling
              Go = single(r_lat_data.Ggrid_RLU(iGo, :));
              phase_shift_coarse = exp(-2 * pi * 1i * R_sampling_RLU * Go.');
              c_rho = c_rho .* phase_shift_coarse;
            elseif use_coarse_G0_phase
              Go = single(r_lat_data.Ggrid_RLU(iGo, :));
              R_coarse_scal = R_coarse_RLU .* single(double(fftgrid_i) ./ double(fftgrid_c));
              phase_shift_coarse = exp(-2 * pi * 1i * R_coarse_scal * Go.');
              c_rho = c_rho .* phase_shift_coarse;
            end

            iscp = [ob, ikpibz, ikprot, ispin];
            param = [];
            param.is = isc;
            param.os = iscp;
            param.qs = [iGo, iqibz, iqrot];
            ngrho_left = SCATTER_Bamp(param);

            Ex_t = sum(vcoul_q .* abs(ngrho_left).^2);
            Ex_ISDF = c_rho' * tildeVq(:, :, iqibz) * c_rho;
            Esum2 = Esum2 + Ex_t^2;
            EsumISDF2 = EsumISDF2 + Ex_ISDF^2;
            DiffEsum2 = DiffEsum2 + (Ex_t - Ex_ISDF)^2;
            Esum = Esum + Ex_t + Ex_ISDF;
            DiffEsum = DiffEsum + abs(Ex_t - Ex_ISDF);
          end
        end % iqbz

        timing.LIVE_timing(1);
      end % ib
    end % ispin
  end % ik
end

function isdftest_ensure_timing_initialized()
  try
    timing.get();
  catch %#ok<CTCH>
    timing.driver();
  end
end

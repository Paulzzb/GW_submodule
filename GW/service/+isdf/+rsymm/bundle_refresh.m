% License-Identifier: GPL
%
% Copyright (C) 2026 ZZ
%
% Last modified: 2026/04/19

function bundle_refresh(id, N_new, indices_new, idnew)
  % Add some new points to the bundle, and refresh the bundle structure.
  %
  %   bundle_refresh(id, N_new, indices_new)
  %   bundle_refresh(id, N_new, indices_new, idnew)
  %
  %   Optional WF consistency checks: isdf.debug.init_from_config (see ISDF_debug.md),
  %   tag rsymm/bundle_refresh.
  %
  % If idnew is omitted or empty, idnew = id (update that slot only).
  % If idnew ~= id, copy the full isdf_m from id to idnew, then refresh
  % bundle_struct on idnew; slot id is left unchanged until overwritten elsewhere.

  if nargin < 1 || isempty(id)
    error('isdf.rsymm:adaptive_isdf_closed_r_bundle:Id', 'ISDF id is required.');
  end

  id = int32(id);
  if nargin < 4 || isempty(idnew)
    idnew = id;
  else
    idnew = int32(idnew);
  end

  rot_tol = 1e-4;

  isdf_data = isdf.get(id);
  if idnew ~= id
    list = isdf.isdf_list();
    if list(idnew).empty
      isdf_dup = isdf_data;
      isdf_dup.id = idnew;
      isdf.save2mod(isdf_dup, idnew);
    end
    isdf_data = isdf.get(idnew);
  end
  isdf_data1 = isdf.get(id);
  fft_data = FFT.get();
  wf_data = wave_functions.get();
  k_data = lattice.manager('k', 'get');
  symm_data = symmetry.get();

  nr = double(fft_data.nr);
  nsym = symm_data.nsym;
  fftgrid = single(fft_data.fftgrid);
  g = double(fftgrid(:)).';
  R_rot = zeros(N_new, nsym);
  


  N_sampling_old = isdf_data.bundle_struct.N_sampling;
  N_sampling_new = N_sampling_old + N_new;
  N_coarse = double(isdf_data.bundle_struct.N_coarse);

  Nb_old = isdf_data.bundle_struct.N_bundle;
  Nb_new = Nb_old;
  Rgrid_b_new = zeros(nr, 3);
  Rgrid_b_new(1:Nb_old, :) = isdf_data.bundle_struct.R_grid_bundle;

  WF_b_new = zeros(size(wf_data.c));
  WF_b_new(1:Nb_old, :, :, :) = isdf_data.bundle_struct.WF_bundle;

  s2b_new = int32(zeros(N_sampling_new, 1));
  s2b_new(1:N_sampling_old) = isdf_data.bundle_struct.sampling2bundle;


  Rrot_b_new = int32(zeros(nr, nsym));
  Rrot_b_new(1:Nb_old, :) = isdf_data.bundle_struct.R_rot_in_bundle;

  new_generator = zeros(N_new, 3);
  N_new_generator = 0;
  r_rot = zeros(nsym, 3);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Expand bundle
  % first, set old bundle --> refreshed_bundle {Rgrid_b_new}
  % For each new point {r_new}, first check if it is in the `refreshed_bundle'
  % if not
  %     use it to calculate all {S*r_new} for all symmetry operations S
  %     generate a new orbit by using `unique' : unique{S*r_new} --> {ind_orbit_finegrid}
  %     add the new orbit to the refreshedbundle :
  %         [Rgrid_b_new; Rgrid_b_new(ind_orbit_finegrid, :)] --> Rgrid_b_new
  % end if
  % add that sampling point indices in the bundle --> s2b_new
  for inew = 1:N_new
    r_new = single( fft_data.Rgrid_RLU(indices_new(inew), :) );
    tmp = find(all( abs( Rgrid_b_new(1:Nb_new, :) - r_new ) < rot_tol, 2));
    if length(tmp) > 2
      error('isdf.rsymm:bundle_refresh:rNewNotInBundle', sprintf( ...
        'R_sampling_RLU row %d is not any row of bundle_struct.R_grid_bundle (single, ==). r_new = [%g %g %g].', ...
        inew, double(r_new(1)), double(r_new(2)), double(r_new(3))));
    end
    if isempty(tmp)
      % If the new point is not in the old bundle, use it to generate a new orbit
      N_new_generator = N_new_generator + 1;
      new_generator(N_new_generator, :) = r_new;
      for isym = 1:nsym
        r_rot(isym, :) = r_new * symm_data.rot_mtrx_RLU_R(:, :, isym);
      end
      r_rot = mod(round(r_rot) + fftgrid, fftgrid);
      ind_rot_finegrid = 1 + r_rot(:, 1) + fftgrid(1) * r_rot(:, 2) + fftgrid(1) * fftgrid(2) * r_rot(:, 3);
      if round(ind_rot_finegrid) - ind_rot_finegrid > rot_tol
        error('isdf.rsymm:bundle_refresh:indRotFinegridNotInteger', sprintf( ...
          'ind_rot_finegrid is not an integer. ind_rot_finegrid = %s.', mat2str(ind_rot_finegrid(:).', 16)));
      end
      ind_rot_finegrid = int32( round(ind_rot_finegrid) );
      ind_orbit_finegrid = unique(ind_rot_finegrid(:), 'stable');
      No = length(ind_orbit_finegrid);
      %
      Rgrid_b_new(Nb_new+1:Nb_new+No, :) = fft_data.Rgrid_RLU(ind_orbit_finegrid, :);
      WF_b_new(Nb_new+1:Nb_new+No, :, :, :) = wf_data.c(ind_orbit_finegrid, :, :, :);
      s2b_new(N_sampling_old + inew) = Nb_new + 1; 

      Nb_new = Nb_new + No;
    else
      s2b_new(N_sampling_old + inew) = tmp(1);
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Calculate the index of the new bundle points in the fine grid
  % If some bundle points are not inside the fine grid, ignore it.
  % ind_b_in_finegrid : Nr x 1, where only 1:Nb_in_finegrid are valid
  %    --> index of the bundle points that are inside the fine grid
  % Nb_in_finegrid    :
  %    --> number of the bundle points that are inside the fine grid
  % ind_b_not_in_finegrid : Nr x 1, where only 1:Nb_not_in_finegrid are valid
  %    --> index of the bundle points that are not inside the fine grid
  % Nb_not_in_finegrid :
  %    --> number of the bundle points that are not inside the fine grid
  % finegrid2newb     : Nr x 1
  %    --> index that map the fine grid index to the bundle index
  % bundle2finegrid   : Nb_new x 1, where only 1:Nb_new are valid
  %    --> index that map the bundle index to the fine grid index
  % NOTE: bundle2finegrid is the inverse of finegrid2newb
  %    --> 
  % 
  ind_b_in_finegrid = zeros(Nb_new, 1);
  ind_b_not_in_finegrid = zeros(Nb_new, 1);
  count1 = 0;
  count2 = 0;
  for irb = 1:Nb_new
    if norm(Rgrid_b_new(irb, :) - round(Rgrid_b_new(irb, :))) < rot_tol
      count1=count1+1;
      ind_b_in_finegrid(count1) = irb;
    else
      count2=count2+1;
      ind_b_not_in_finegrid(count2) = irb;
    end
  end
  Nb_in_finegrid = count1;
  Nb_not_in_finegrid = count2;


  finegrid2newb = zeros(nr, 1);
  irbninfft_i = 1 + Rgrid_b_new(ind_b_in_finegrid(1:Nb_in_finegrid), 1) ...
               + fftgrid(1) * Rgrid_b_new(ind_b_in_finegrid(1:Nb_in_finegrid), 2) ...
               + fftgrid(1) * fftgrid(2) * Rgrid_b_new(ind_b_in_finegrid(1:Nb_in_finegrid), 3);
  finegrid2newb(irbninfft_i) = ind_b_in_finegrid(1:Nb_in_finegrid);

  bundle2finegrid = int32(round( ...
    1 + Rgrid_b_new(ind_b_in_finegrid(1:Nb_in_finegrid), 1) ...
    + fftgrid(1) * Rgrid_b_new(ind_b_in_finegrid(1:Nb_in_finegrid), 2) ...
    + fftgrid(1) * fftgrid(2) * Rgrid_b_new(ind_b_in_finegrid(1:Nb_in_finegrid), 3)));

  % Integer-RLU bundle rows must match IBZ grid samples wf_data.c(lin,:,:,:).
  % Old WF_bundle rows often come from coarse / ISDF coefficient paths and can
  % disagree with WF_apply_symm by ~1e-2--1e-3 here, which poisons downstream checks.
  % for ksync = 1:Nb_in_finegrid
  %   irb = double(ind_b_in_finegrid(ksync));
  %   lin = double(bundle2finegrid(ksync));
  %   WF_b_new(irb, :, :, :) = wf_data.c(lin, :, :, :);
  % end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Calculate the rotation in bundle
  % Rrot_b_new     : Nb_new x nsym
  % Rgrid_b_new(Rrot_b_new(irb, isym), :) = Rgrid_b_new(irb, :) * rot_mtrx(:, :, isym)
  %
  % First for `bundle points' that are inside the fine grid
  Rgrid_b_not_in_finegrid = Rgrid_b_new(ind_b_not_in_finegrid(1:Nb_not_in_finegrid), :);
  for isym = 1:nsym
    rot_mtrx = symm_data.rot_mtrx_RLU_R(:, :, isym);
    % Integer lattice: R*S should land on Z^3 before mod (same idea as FFT/driver.m).
    % Use double(R) to avoid false rejects from single*double roundoff; tolerance 1e-3 matches FFT/driver.
    M2_r_RLU = double(Rgrid_b_new(ind_b_in_finegrid(1:Nb_in_finegrid), :)) * rot_mtrx;
    res_int = abs(M2_r_RLU - round(M2_r_RLU));
    if max(res_int(:)) > rot_tol
      error('isdf.rsymm:bundle_refresh:RotInBundleNotInteger', sprintf( ...
        'M2_r_RLU is not near-integer (max|R*S-round(R*S)|=%g). M2_r_RLU = %s.', ...
        max(res_int(:)), mat2str(M2_r_RLU(:).', 16)));
    end
    M2_r_RLU = round(M2_r_RLU);
    iv_mod = double(mod(M2_r_RLU + fftgrid, fftgrid));
    indR = 1 + iv_mod(:,1) + iv_mod(:,2)*fftgrid(1) ...
             + iv_mod(:,3)*fftgrid(1)*fftgrid(2);
    Rrot_b_new(ind_b_in_finegrid(1:Nb_in_finegrid), isym) = finegrid2newb(indR);
  end
  % Then for `bundle points' that are not inside the fine grid
  % Notice: no integer lattice guarantee; debug-print rows where mod(residual,G) is ~0 in all 3 comps.
  if Nb_not_in_finegrid > 0
    for isym = 1:nsym
      rot_mtrx = symm_data.rot_mtrx_RLU_R(:, :, isym);
      M2_r_RLU = double(Rgrid_b_new(ind_b_not_in_finegrid(1:Nb_not_in_finegrid), :)) * rot_mtrx;
      diff_mat = M2_r_RLU - double(Rgrid_b_not_in_finegrid);
      tmp = mod(diff_mat + double(fftgrid), double(fftgrid));
      tol_zero = rot_tol;
      row_all_near_zero = all(abs(tmp) < tol_zero, 2);
      iz = find(row_all_near_zero);
      if ~isempty(iz)
        irb_bundle = ind_b_not_in_finegrid(iz);
        fprintf(1, ['[isdf.rsymm.bundle_refresh] isym=%d: tmp rows with all 3 comps |.|<%g ', ...
          '(subset ix / bundle irb):\n  ix=%s\n  irb=%s\n'], ...
          isym, tol_zero, mat2str(iz(:).'), mat2str(irb_bundle(:).'));
      end
    end
  end

  if isdf.debug.on('rsymm/bundle_refresh')
    % --- verify the rotation (torus consistency) ---
    for isym = 1:nsym
      rot_mtrx = symm_data.rot_mtrx_RLU_R(:, :, isym);
      for irb = ind_b_in_finegrid(1:Nb_in_finegrid)
        rbg = Rgrid_b_new(irb, :);
        rbg_rot = mod(rbg * rot_mtrx + fftgrid, fftgrid);
        irbrot = Rrot_b_new(irb, isym);
        rbg_rot_tgt = single(Rgrid_b_new(irbrot, :));
        dmx = max(abs(rbg_rot - rbg_rot_tgt));
        msg = sprintf(['Rotation map mismatch: irb=%d isym=%d irbrot=%d; ', ...
          'max|mod(R*S,G) - R(irbrot,:)| = %g.'], irb, isym, irbrot, dmx);
        isdf.debug.react(dmx > rot_tol, msg, 'bundle_refresh_rot_v0');
      end
    end

    % --- discrete rotation map on bundle rows (integer RLU) ---
    gmod = double(fftgrid);
    for isym = 1:nsym
      rot_mtrx = symm_data.rot_mtrx_RLU_R(:, :, isym);
      for k = 1:Nb_in_finegrid
        irb = double(ind_b_in_finegrid(k));
        irbrot = double(Rrot_b_new(irb, isym));
        bad_idx = irbrot < 1 || irbrot > Nb_new || irbrot ~= floor(irbrot);
        msg_idx = sprintf('R_rot_in_bundle(%d,%d)=%g is not an integer in [1,%d].', ...
          irb, isym, irbrot, Nb_new);
        isdf.debug.react(bad_idx, msg_idx, 'bundle_refresh_rot_idx');
        if bad_idx
          continue
        end
        v_rot = mod(round(double(Rgrid_b_new(irb, :)) * rot_mtrx) + gmod, gmod);
        v_tgt = double(Rgrid_b_new(irbrot, :));
        dgeom = max(abs(v_rot - v_tgt));
        msg_geom = sprintf(['Rotation map mismatch: irb=%d isym=%d irbrot=%d; ', ...
          'max|mod(round(R*S),G) - R(irbrot,:)| = %g.'], irb, isym, irbrot, dgeom);
        isdf.debug.react(dgeom > rot_tol, msg_geom, 'bundle_refresh_rot_geom');
      end
    end

    nspin = wf_data.nspin;
    nbz = k_data.nbz;
    is_t_rev = symm_data.is_t_rev;
    R_sampling = Rgrid_b_new(s2b_new(N_coarse+1:N_sampling_new), :);

    for ispin = 1:nspin
      for ikbz = 1:nbz
        ikibz = k_data.bz2ibz(ikbz);
        ikrot = k_data.bz2rot(ikbz);
        for ib = 1:nb
          isc = [ib, ikibz, ikrot, ispin];
          wf = wave_functions.WF_apply_symm(isc);
          wf_in_bundle = wf(bundle2finegrid);
          wf_ikibz_b = WF_b_new(:, ib, ikibz, ispin);
          if ikrot > nsym / (is_t_rev + 1)
            ikrot_wf = ikrot - nsym / (is_t_rev + 1);
            inv_ikrot_wf = symm_data.inv_rot_index(ikrot_wf);
            conjflag = true;
          else
            ikrot_wf = ikrot;
            inv_ikrot_wf = symm_data.inv_rot_index(ikrot_wf);
            conjflag = false;
          end
          indinbundle = Rrot_b_new(1:Nb_new, inv_ikrot_wf);
          wf_ikbz_b = wf_ikibz_b(indinbundle);
          if conjflag
            wf_ikbz_b = conj(wf_ikbz_b);
          end
          isdf.debug.react(norm(wf_in_bundle - wf_ikbz_b(ind_b_in_finegrid(1:Nb_in_finegrid))) > 8e-5, ...
            'The wavefunction on the bundle is not correct', 'bundle_refresh_wf_bundle');
          wf_sampling1 = wf(indices_new(1:N_new));
          wf_sampling2 = zeros(N_sampling_new - N_coarse, 1);
          Srotk_R_sampling = R_sampling * symm_data.rot_mtrx_RLU_R(:, :, inv_ikrot_wf);
          Srotk_R_sampling = mod(round(Srotk_R_sampling) + fftgrid, fftgrid);
          ind_Srotk_R_sampling = 1 + Srotk_R_sampling(:, 1) + fftgrid(1) * Srotk_R_sampling(:, 2) ...
                                 + fftgrid(1) * fftgrid(2) * Srotk_R_sampling(:, 3);
          ind_Srotk_R_sampling = int32(round(ind_Srotk_R_sampling));
          ind_in_bundle = finegrid2newb(ind_Srotk_R_sampling);
          wf_sampling2 = wf_ikibz_b(ind_in_bundle);
          if conjflag
            wf_sampling2 = conj(wf_sampling2);
          end
          wf_sampling3 = wf_ikbz_b(s2b_new(N_coarse+1:N_sampling_new));
          isdf.debug.react(norm(wf_sampling1 - wf_sampling2) > 8e-5, ...
            'The wavefunction on the sampling is not correct', 'bundle_refresh_sampling12');
          isdf.debug.react(norm(wf_sampling1 - wf_sampling3) > 8e-5, ...
            'The wavefunction on the sampling is not correct', 'bundle_refresh_sampling13');
        end
      end
    end

    for ispin = 1:nspin
      for ikbz = 1:nbz
        ikibz = k_data.bz2ibz(ikbz);
        ikrot = k_data.bz2rot(ikbz);
        if ikrot > nsym / (is_t_rev + 1)
          ikrot_wf = ikrot - nsym / (is_t_rev + 1);
          inv_ikrot_wf = symm_data.inv_rot_index(ikrot_wf);
        else
          ikrot_wf = ikrot;
          inv_ikrot_wf = symm_data.inv_rot_index(ikrot_wf);
        end
        for ib = 3:5
          isc = [ib, ikibz, ikrot, ispin];
          wf = wave_functions.WF_apply_symm(isc);
          indinbundle = Rrot_b_new(1:Nb_new, inv_ikrot_wf);
          ind_sampling_fine = s2b_new(N_coarse+1:N_sampling_new);
          wf_sampling = wf(indinbundle(ind_sampling_fine));

          R_sampling_b = Rgrid_b_new(s2b_new(N_coarse+1:N_sampling_new), :);
          Srotk_R_sampling = R_sampling_b * symm_data.rot_mtrx_RLU_R(:, :, inv_ikrot_wf);
          Srotk_R_sampling = mod(round(Srotk_R_sampling) + fftgrid, fftgrid);
          ind_Srotk_R_sampling = 1 + Srotk_R_sampling(:, 1) + fftgrid(1) * Srotk_R_sampling(:, 2) ...
                                 + fftgrid(1) * fftgrid(2) * Srotk_R_sampling(:, 3);
          ind_Srotk_R_sampling = int32(round(ind_Srotk_R_sampling));
          ind_bundle = finegrid2newb(ind_Srotk_R_sampling);
          wf_b = WF_b_new(:, ib, ikibz, ispin);
          wf_Srotk_R_sampling = wf_b(ind_bundle);
          isdf.debug.react(norm(wf_sampling - wf_Srotk_R_sampling) > 8e-5, ...
            'The wavefunction on the bundle is not correct (sampling vs rotated path, ib=3:5).', ...
            'bundle_refresh_wf_sampling_rotate');
        end
      end
    end
  end



  bs_new = struct();
  bs_new.N_bundle = Nb_new;
  bs_new.N_sampling = N_sampling_new;
  bs_new.R_grid_bundle = Rgrid_b_new(1:Nb_new, :);
  bs_new.sampling2bundle = s2b_new;
  bs_new.R_rot_in_bundle = Rrot_b_new(1:Nb_new, :);
  bs_new.WF_bundle = WF_b_new(1:Nb_new, :, :, :);
  bs_new.N_coarse = isdf_data.bundle_struct.N_coarse;
  isdf_data.bundle_struct = bs_new;
  %
  isdf.save2mod(isdf_data, idnew);

  if isdf.debug.on('rsymm/bundle_refresh')
    % Post-save bundle test: WF on bs_new vs WF_apply_symm; S_q on bundle vs direct.
    is_t_rev = symm_data.is_t_rev;
    % Same integer-RLU bundle rows as finegrid2newb / first verification block
    bundle2finegrid = int32(round( ...
      1 + bs_new.R_grid_bundle(ind_b_in_finegrid(1:Nb_in_finegrid), 1) ...
      + fftgrid(1) * bs_new.R_grid_bundle(ind_b_in_finegrid(1:Nb_in_finegrid), 2) ...
      + fftgrid(1) * fftgrid(2) * bs_new.R_grid_bundle(ind_b_in_finegrid(1:Nb_in_finegrid), 3)));
    for ispin = 1:wf_data.nspin
      for ikbz = 1:k_data.nbz
        for ib = 3:5
          ikibz = k_data.bz2ibz(ikbz);
          ikrot = double(k_data.bz2rot(ikbz));
          wf_ibz_in_bundle_a = bs_new.WF_bundle(:, ib, ikibz, ispin);
          if ikrot > nsym / (is_t_rev + 1)
            ikrot_wf = ikrot - nsym / (is_t_rev + 1);
            ikrot_wf = symm_data.inv_rot_index(ikrot_wf);
            conjflag = true;
          else
            ikrot_wf = ikrot;
            ikrot_wf = symm_data.inv_rot_index(ikrot_wf);
            conjflag = false;
          end
          indices = bs_new.R_rot_in_bundle(:, ikrot_wf);
          wf_bz_in_bundle = wf_ibz_in_bundle_a(indices);
          if conjflag
            wf_bz_in_bundle = conj(wf_bz_in_bundle);
          end
          wf_tmp = wf_bz_in_bundle(1:Nb_new);
          wf_tmp_finegrid = wf_tmp(ind_b_in_finegrid(1:Nb_in_finegrid));
          isc = [ib, ikibz, ikrot, ispin];
          wf_dir = wave_functions.WF_apply_symm(isc);
          wf_bz_in_bundle_dir = wf_dir(bundle2finegrid);
          nmf = norm(wf_tmp_finegrid - wf_bz_in_bundle_dir);
          msg_wf = sprintf(['WF on bundle vs WF_apply_symm disagree. isc=[%d %d %d %d] ', ...
            'norm(fine diff)=%g norm(full bundle vs dir)=%g.'], ib, ikibz, ikrot, ispin, nmf, ...
            norm(wf_bz_in_bundle - wf_bz_in_bundle_dir));
          isdf.debug.react(nmf > 8e-5, msg_wf, 'bundle_refresh_postsave_wf_dir');
          for iq = 1:k_data.nbz
            iqrot = double(k_data.bz2rot(iq));
            ind = bs_new.R_rot_in_bundle(:, iqrot);
            wf_Sq_bundle = wf_bz_in_bundle(ind);
            wf_Sq_xalpha = wf_Sq_bundle(bs_new.sampling2bundle(1:N_sampling_new));
            wf_Sq_xalpha_fine = wf_Sq_xalpha(N_coarse+1:N_sampling_new);
            wf_Sq_xalpha_coarse = wf_Sq_xalpha(1:N_coarse);
            R_grid_sampling = bs_new.R_grid_bundle(bs_new.sampling2bundle(N_coarse+1:N_sampling_new), :);
            Sq_R_grid_sampling = R_grid_sampling * symm_data.rot_mtrx_RLU_R(:, :, iqrot);
            nrint = norm(Sq_R_grid_sampling - round(Sq_R_grid_sampling));
            msg_rot = sprintf( ...
              'Sampling rotation on bundle not near-integer grid (ib=%d ikibz=%d iq=%d ispin=%d nrint=%g).', ...
              ib, ikibz, iq, ispin, nrint);
            isdf.debug.react(nrint > 5e-4, msg_rot, 'bundle_refresh_postsave_sampling_rot');
            Sq_R_grid_sampling = single(round(Sq_R_grid_sampling));
            Sq_R_grid_sampling = mod(round(Sq_R_grid_sampling) + fftgrid, fftgrid);
            ind = 1 + Sq_R_grid_sampling(:, 1) + fftgrid(1) * Sq_R_grid_sampling(:, 2) ...
                   + fftgrid(1) * fftgrid(2) * Sq_R_grid_sampling(:, 3);
            wf_Sq_xalpha_fine_dir = wf_dir(ind);
            nmfine = norm(wf_Sq_xalpha_fine - wf_Sq_xalpha_fine_dir);
            msg_fine = sprintf( ...
              'Fine-grid WF after S_q disagree (ib=%d ikibz=%d iq=%d ispin=%d norm=%g).', ...
              ib, ikibz, iq, ispin, nmfine);
            isdf.debug.react(nmfine > 8e-5, msg_fine, 'bundle_refresh_postsave_Sq_fine');
            wf_ikibz_coarse = isdf_data1.coeff_seper(:, ib, ikibz, ispin);
            wf_ikibz_coarse = wf_ikibz_coarse(1:N_coarse);
            wf_ikbz_coarse = isdf.coeff.isdf_apply_symm_on_coarse(id, wf_ikibz_coarse, ikrot);
            ind_Sq_R_coarse = isdf_data1.R_rot_coarse(:, iqrot);
            wf_Sq_xalpha_coarse_dir = wf_ikbz_coarse(ind_Sq_R_coarse);
            nmcoarse = norm(wf_Sq_xalpha_coarse - wf_Sq_xalpha_coarse_dir);
            msg_coarse = sprintf( ...
              'Coarse WF after S_q disagree (ib=%d ikibz=%d iq=%d ispin=%d norm=%g).', ...
              ib, ikibz, iq, ispin, nmcoarse);
            isdf.debug.react(nmcoarse > 8e-5, msg_coarse, 'bundle_refresh_postsave_Sq_coarse');
          end
        end
      end
    end
  end % if debug

end

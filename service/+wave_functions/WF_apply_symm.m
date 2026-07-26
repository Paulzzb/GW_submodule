function WF_symm = WF_apply_symm(isc)
  % ELsh, spin_sop, idt_index, myid, 
  persistent firsttime nsym is_t_rev fftgrid nr

  if nargin == 1 && (ischar(isc) || (isstring(isc) && isscalar(isc)))
    cmd = lower(string(isc));
    if cmd == "reset"
      firsttime = [];
      nsym = [];
      is_t_rev = [];
      fftgrid = [];
      nr = [];
      WF_symm = [];
      return
    end
  end

  if isempty(firsttime)
    firsttime = true;
  end

  if firsttime
    symm_data = symmetry.manager('get');
    nsym = symm_data.nsym;
    is_t_rev= symm_data.is_t_rev;
    fft_data = FFT.manager('get');
    fftgrid = fft_data.fftgrid;
    nr = prod(fftgrid);
    %
    firsttime = false;
  end

  WF = wave_functions.manager('get');
  
  ib = isc(1);
  ikibz = isc(2);
  ispin = isc(4);
  isymm = isc(3);

  WF_symm = complex(zeros(nr, WF.n_spinor, 'double'));

  if isymm == 1
    WF_symm = WF.c(:, ib, ikibz, ispin);
    return
  end

  fft_data = FFT.manager('get');
  symm_data = symmetry.manager('get');

  if WF.n_spinor == 1
    inv_isymm = symm_data.inv_rot_index(isymm);
    ind = fft_data.R_rot(:, inv_isymm);
    WF_symm(:, 1) = WF.c(ind, ib, ikibz, ispin);

  elseif WF.n_spinor == 2
    error("Support for spinor wave function is developing");
    % ind = fft_data.R_rot(:, isymm);

    % WF_symm(:, 1) = spin_sop(1, 1, isymm) .* WF.c(ind, ib, ikibz, ispin) + ...
    %                 spin_sop(1, 2, isymm) .* WF.c(ind, ib, ikibz, ispin);

    % WF_symm(:, 2) = spin_sop(2, 1, isymm) .* WF.c(ind, ib, ikibz, ispin) + ...
    %                 spin_sop(2, 2, isymm) .* WF.c(ind, ib, ikibz, ispin);
  else
    error('WF_apply_symm:UnsupportedSpinor', ...
      'Unsupported n_spinor = %d', WF.n_spinor);
  end

  if isymm > nsym / (is_t_rev + 1)
    WF_symm = conj(WF_symm);
  end
end

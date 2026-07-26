% License-Identifier: GPL
%
% Copyright (C) 2026 ZZ
%
% Last modified: 2026/04/19

function out = adaptive_isdf_closed_r_bundle(id)
%ADAPTIVE_ISDF_CLOSED_R_BUNDLE  Adaptive ISDF: symmetry-closed r-grid subset, wf on IBZ, R_rot on subset.
%
%   out = isdftest.rsymm.adaptive_isdf_closed_r_bundle(id)
%
% For interp_scheme == "adaptive" only. Builds the smallest subset U of fine FFT
% sites (linear indices 1..nr) that contains every ISDF centroid and is closed
% under all columns of FFT.R_rot (i.e. r' = S r on the discrete torus), by BFS
% on the same graph as adaptiveisdf's orbit report.
%
% Wavefunctions are taken only on the IBZ mesh (indices 1..nibz), as stored in
% wave_functions.c(:, ib, ik_ibz, ispin).
%
% Outputs (struct):
%   lin_fine     鈥?N_sub x 1 int32, sorted unique fine-grid linear indices in U
%   wf           鈥?complex double, size [N_sub, nb, nibz, nspin]
%   R_rot_sub    鈥?int32 [N_sub, nsym], subset-index rotation (see comments near out.*)
%   n_sub, nr, nsym 鈥?scalars (|U|, nr, number of symmetry columns in R_rot)
%   centroid_lin 鈥?fine indices for rows in Nrange (see body: adaptive uses extra rows only)
%
% Notes:
%   - BFS seeds are unique(indices_generater) where indices_generater maps R_sampling_RLU(Nrange,:)
%     to fine-grid lines (non-adaptive: all rows; adaptive: N_coarse+1 : N_coarse+N_extra).

  if nargin < 1 || isempty(id)
    error('isdftest.rsymm:adaptive_isdf_closed_r_bundle:Id', 'ISDF id is required.');
  end

  isdf_data = isdftest.get(id);
  fft_data = FFT.get();
  wf_data = wave_functions.get();
  k_data = lattice.manager('k', 'get');

  R_rot = fft_data.R_rot;
  nr = double(fft_data.nr);
  nsym = int32(size(R_rot, 2));
  if nsym < 1
    error('isdftest.rsymm:adaptive_isdf_closed_r_bundle:R_rot', ...
    'FFT.R_rot has no symmetry columns.');
  end

  Nisdf = double(isdf_data.nisdf);
  if Nisdf < 1
    error('isdftest.rsymm:adaptive_isdf_closed_r_bundle:Empty',...
         'Empty ISDF sampling (nisdf<1).');
  end

  if isdf_data.interp_scheme ~= "adaptive"
    N_old = 0;
    N_new = Nisdf;
  else
    N_old = isdf_data.N_coarse;
    N_new = isdf_data.N_extra;
  end
  Nrange = N_old+1:N_new+N_old;

  indices_generater = r_sampling_rows_to_lin(isdf_data, fft_data, Nrange);
  seeds = unique(indices_generater(:), 'stable');

  U = orbit_union_bfs_from_seeds(seeds, double(R_rot), nr);
  lin_fine = int32(find(U));
  lin_fine = sort(lin_fine(:));
  n_sub = int32(numel(lin_fine));

  full_to_sub = zeros(nr, 1, 'int32');
  full_to_sub(lin_fine) = int32(1:double(n_sub));

  R_rot_sub = zeros(double(n_sub), double(nsym), 'int32');
  for is = 1:double(nsym)
    dest = full_to_sub(R_rot(lin_fine, is));
    if any(dest == 0)
      error('isdftest.rsymm:adaptive_isdf_closed_r_bundle:NotClosed', ...
        'Orbit union not closed under symmetry %d (check FFT.R_rot).', is);
    end
    R_rot_sub(:, is) = dest;
  end

  nb = int32(wf_data.nb);
  nibz = int32(k_data.nibz);
  nspin = int32(wf_data.nspin);

  wf = wf_data.c(lin_fine, 1:nb, 1:nibz, 1:nspin);

  % --- Pack outputs (see header block for formulas) ---
  % lin_fine: 闂寘 U 鍐呮墍鏈夌粏缃戞牸鐐圭殑绾挎寚鏍?1..nr锛屽崌搴忋€佸幓閲嶏紱闀垮害 n_sub銆?  out.lin_fine = lin_fine;
  % wf: 娉㈠嚱鏁颁粎鍦ㄩ棴鍖呯偣 lin_fine 涓婂彇鍊硷紱缁存暟 [n_sub, nb, nibz, nspin]锛屼笌 wave_functions.c 鐨?IBZ 缁翠竴鑷淬€?  out.wf = wf;
  % R_rot_sub(i_sub,is): 瀛愰泦缂栧彿涓嬶紝绗?i_sub 涓棴鍖呯偣鍦ㄧ┖闂村绉?is 涓嬬殑鍍?鈫?瀛愰泦缂栧彿 j_sub锛?  %   婊¤冻 lin_fine(j_sub) = R_rot(lin_fine(i_sub), is)銆?  out.R_rot_sub = R_rot_sub;
  % n_sub: 闂寘鐐逛釜鏁?|U|銆?  out.n_sub = n_sub;
  % nr: 缁?FFT 瀹炵┖闂存€绘牸鐐规暟锛坧rod(fftgrid)锛夈€?  out.nr = int32(nr);
  % nsym: FFT.R_rot 鐨勫垪鏁帮紙浠ｇ爜閲岀敤鍒扮殑绌洪棿瀵圭О鎿嶄綔涓暟锛夈€?  out.nsym = nsym;
  % centroid_lin: 褰撳墠 Nrange 瀵瑰簲閭ｄ簺 R_sampling 琛屾槧鍒扮粏缃戞牸鐨勭嚎鎸囨爣锛堜笌涓婇潰 seeds 鍚屾簮锛涜搴忓悓 coeff_seper 涓娈碉級銆?  out.centroid_lin = indices_generater;
end

% --- Same mapping as adaptiveisdf_r_sampling_rows_to_lin (adaptiveisdf.m) ---

function lin = r_sampling_rows_to_lin(isdf_data, fft_data, Nrange)
  ni = double(fft_data.fftgrid(:)).';
  Rgrid = double(fft_data.Rgrid_RLU);
  Nmu = length(Nrange);
  lin = zeros(Nmu, 1);
  Rs = double(isdf_data.R_sampling_RLU(Nrange, :));
  for i = 1:Nmu
    v = round(Rs(i, :));
    v = mod(v, ni);
    [is_hit, k] = ismember(v, Rgrid, 'rows');
    if ~is_hit
      dd = zeros(size(Rgrid, 1), 3);
      for ddim = 1:3
        t = abs(Rgrid(:, ddim) - v(ddim));
        dd(:, ddim) = min(t, min(abs(t - ni(ddim)), abs(t + ni(ddim))));
      end
      [~, k] = min(sum(dd, 2));
    end
    lin(i) = k;
  end
end

function U = orbit_union_bfs_from_seeds(seeds, R_rot, nr)
  U = false(nr, 1);
  for k = 1:numel(seeds)
    U = U | orbit_mask_bfs(seeds(k), R_rot, nr);
  end
end

function Om = orbit_mask_bfs(seed, R_rot, nr)
  Om = false(nr, 1);
  dq = seed;
  Om(seed) = true;
  head = 1;
  while head <= numel(dq)
    i = dq(head);
    head = head + 1;
    for is = 1:size(R_rot, 2)
      j = R_rot(i, is);
      if j >= 1 && j <= nr && ~Om(j)
        Om(j) = true;
        dq(end + 1) = j;
      end
    end
  end
end

% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/16

function isdf_vitualize()
% Visualize ISDF interpolation sites in the primitive cell (no return value, no console output).
% Data are read only via service getters: FFT.get, lattice.manager, isdf.get.

  fft_data = FFT.get();
  d_lat = lattice.manager('d_lat', 'get');
  isdf_data = isdf.get();

  A = double(d_lat.a1a2a3);
  fftgrid = double(fft_data.fftgrid(:)).';

  fig = figure('Color', 'w', 'Name', 'ISDF visualize');
  ax = axes(fig);
  hold(ax, 'on');
  axis(ax, 'equal');
  grid(ax, 'on');
  xlabel(ax, 'Cartesian (same units as lattice vectors)');
  ylabel(ax, 'y');
  zlabel(ax, 'z');
  title(ax, sprintf('ISDF id=%d  nisdf=%d  scheme=%s', ...
    double(isdf_data.id), double(isdf_data.nisdf), char(string(isdf_data.interp_scheme))));

  local_draw_unit_cell(ax, A);

  ap = d_lat.atom_pos;
  if ~isempty(ap)
    nat = size(ap, 1);
    xyz_a = reshape(double(ap), nat, 3);
    syms = d_lat.atom_symbol;
    if numel(syms) ~= nat
      syms = repmat({''}, nat, 1);
    end
    [uSym, ~, ic] = unique(string(syms), 'stable');
    cmap = lines(max(7, numel(uSym)));
    for iu = 1:numel(uSym)
      m = ic == iu;
      lab = strtrim(char(uSym(iu)));
      if isempty(lab)
        lab = '?';
      end
      scatter3(ax, xyz_a(m, 1), xyz_a(m, 2), xyz_a(m, 3), 120, cmap(iu, :), 'o', ...
        'filled', 'MarkerFaceAlpha', 0.45, 'MarkerEdgeColor', cmap(iu, :) * 0.6, ...
        'DisplayName', sprintf('atoms: %s', lab));
    end
  end

  R = isdf_data.R_sampling_RLU;
  nisdf = double(isdf_data.nisdf);
  if isempty(R) || nisdf < 1
    legend(ax, 'show', 'Location', 'bestoutside');
    view(ax, 35, 22);
    return
  end

  ns = min(nisdf, size(R, 1));
  R = R(1:ns, :);
  cart_s = local_rlu_rows_to_cart(R, fftgrid, A);

  scheme = lower(strtrim(char(string(isdf_data.interp_scheme))));
  Nco = double(isdf_data.N_coarse);
  Nex = double(isdf_data.N_extra);

  if strcmp(scheme, 'adaptive') && Nex > 0 && Nco > 0
    n_co = min(Nco, ns);
    n_ex_show = min(max(0, ns - n_co), Nex);
    scatter3(ax, cart_s(1:n_co, 1), cart_s(1:n_co, 2), cart_s(1:n_co, 3), 36, [0.2 0.55 0.95], ...
      'filled', 'MarkerFaceAlpha', 0.85, 'DisplayName', 'isdf coarse');
    if n_ex_show > 0
      j1 = n_co + 1;
      j2 = n_co + n_ex_show;
      scatter3(ax, cart_s(j1:j2, 1), cart_s(j1:j2, 2), cart_s(j1:j2, 3), 36, [0.95 0.35 0.1], ...
        'filled', 'MarkerFaceAlpha', 0.9, 'DisplayName', 'isdf added');
    end
  elseif strcmp(scheme, 'adaptive') && Nex > 0 && Nco == 0
    scatter3(ax, cart_s(:, 1), cart_s(:, 2), cart_s(:, 3), 36, [0.95 0.35 0.1], ...
      'filled', 'MarkerFaceAlpha', 0.9, 'DisplayName', 'isdf added');
  elseif strcmp(scheme, 'adaptive') && Nco > 0
    scatter3(ax, cart_s(:, 1), cart_s(:, 2), cart_s(:, 3), 36, [0.2 0.55 0.95], ...
      'filled', 'MarkerFaceAlpha', 0.85, 'DisplayName', 'isdf coarse');
  else
    scatter3(ax, cart_s(:, 1), cart_s(:, 2), cart_s(:, 3), 36, [0.15 0.75 0.35], ...
      'filled', 'MarkerFaceAlpha', 0.85, 'DisplayName', 'isdf coarse grid');
  end

  legend(ax, 'show', 'Location', 'bestoutside');
  view(ax, 35, 22);
end

function cart = local_rlu_rows_to_cart(R_rlu, fftgrid, A)
  frac = double(R_rlu) ./ fftgrid;
  cart = frac * A';
end

function local_draw_unit_cell(ax, A)
  c1 = A(:, 1)';
  c2 = A(:, 2)';
  c3 = A(:, 3)';
  V = [0, 0, 0; c1; c2; c3; c1 + c2; c1 + c3; c2 + c3; c1 + c2 + c3];
  E = [1, 2; 1, 3; 1, 4; 2, 5; 2, 6; 3, 5; 3, 7; 4, 6; 4, 7; 5, 8; 6, 8; 7, 8];
  col = [0.1 0.1 0.1];
  for ie = 1:size(E, 1)
    p = V(E(ie, :), :);
    plot3(ax, p(:, 1), p(:, 2), p(:, 3), '-', 'Color', col, 'LineWidth', 1.4, ...
      'HandleVisibility', 'off');
  end
  plot3(ax, nan, nan, nan, '-', 'Color', col, 'LineWidth', 1.4, 'DisplayName', 'primitive cell');
end

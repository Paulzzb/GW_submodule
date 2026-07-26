function E = report_ex_benchmark(E, case_dir, tag)
%REPORT_EX_BENCHMARK  Print/write Ex tables, config.ISDF, and save plots.
%
%   Writes case_dir/Ex.dat (same layout as console report).
%   config.ISDF: inv_strategy and order listed first, then remaining fields.

  if nargin < 3 || isempty(tag)
    tag = local_case_tag(case_dir);
  end
  if nargin < 2 || isempty(case_dir)
    case_dir = pwd;
  end

  E.case_dir = case_dir;
  E.case_tag = tag;
  E = local_write_ex_report(E, case_dir, tag);
  E = local_plot_ex_benchmark(E, case_dir, tag);
end

function E = local_write_ex_report(E, case_dir, tag)
  E.ex_dat = fullfile(case_dir, 'Ex.dat');

  if ~isfield(E, 'Ex_dir') || isempty(E.Ex_dir)
    fprintf('\n[%s] Ex array is empty; skip report.\n', tag);
    return
  end

  fid = fopen(E.ex_dat, 'w');
  if fid < 0
    error('report_ex_benchmark:ExDat', 'Cannot open %s for writing.', E.ex_dat);
  end
  cleanup = onCleanup(@() fclose(fid));

  emit(fid, '[%s] Ex benchmark report', tag);
  emit(fid, 'Generated: %s', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
  emit(fid, 'case_dir: %s', case_dir);
  emit(fid, '');

  if isfield(E, 'config_isdf') && isstruct(E.config_isdf) ...
      && isfield(E.config_isdf, 'ISDF')
    local_emit_isdf_config(fid, E.config_isdf.ISDF);
    emit(fid, '');
  end

  Ex_dir = E.Ex_dir;
  Ex_isdf = E.Ex_isdf;
  dEx = E.dEx;
  [nb, nk] = size(Ex_dir);
  ib = local_band_indices(E, nb);

  emit(fid, '[%s] Ex comparison by ik (eV), size [%d x %d] (ib x ik)', tag, nb, nk);
  emit(fid, '  ||Ex_dir||_F  = %.6e', norm(Ex_dir(:), 'fro'));
  emit(fid, '  ||Ex_isdf||_F = %.6e', norm(Ex_isdf(:), 'fro'));
  emit(fid, '  ||dEx||_F     = %.6e', norm(dEx(:), 'fro'));
  if norm(Ex_dir(:), 'fro') > 0
    emit(fid, '  rel||dEx||_F  = %.6e', norm(dEx(:), 'fro') / norm(Ex_dir(:), 'fro'));
  end

  if isfield(E, 'wall_input')
    emit(fid, '  wall_input      = %.3f s', E.wall_input);
  end
  if isfield(E, 'wall_ex_isdf')
    emit(fid, '  wall_ex_isdf    = %.3f s', E.wall_ex_isdf);
  end
  if isfield(E, 'wall_ex_dir')
    emit(fid, '  wall_ex_dir     = %.3f s', E.wall_ex_dir);
  end

  for ik = 1:nk
    emit(fid, '');
    emit(fid, 'ik = %d', ik);
    emit(fid, '  %4s %12s %12s %12s', 'nb', 'Ex_dir', 'Ex_isdf', 'dEx');
    for ib_out = 1:nb
      emit(fid, '  %4d %12.6f %12.6f %12.6f', ib(ib_out), ...
        Ex_dir(ib_out, ik), Ex_isdf(ib_out, ik), dEx(ib_out, ik));
    end
    emit(fid, '***');
  end

  fprintf('\n[%s] Wrote report: %s\n', tag, E.ex_dat);
end

function local_emit_isdf_config(fid, isdf_cfg)
  emit(fid, '----------- config.ISDF -----------');
  priority = {'inv_strategy', 'order'};
  for k = 1:numel(priority)
    fn = priority{k};
    if isfield(isdf_cfg, fn)
      emit(fid, '%s', local_fmt_isdf_field(fn, isdf_cfg.(fn)));
    end
  end
  names = fieldnames(isdf_cfg);
  for i = 1:numel(names)
    if ismember(names{i}, priority)
      continue
    end
    emit(fid, '%s', local_fmt_isdf_field(names{i}, isdf_cfg.(names{i})));
  end
end

function line = local_fmt_isdf_field(name, val)
  label = sprintf(' %-22s : ', name);
  if islogical(val)
    if isscalar(val)
      line = [label, sprintf('%d', val)];
    else
      line = [label, mat2str(val)];
    end
  elseif isnumeric(val)
    if isscalar(val)
      if mod(val, 1) == 0
        line = [label, sprintf('%d', val)];
      else
        line = [label, sprintf('%.6g', val)];
      end
    else
      line = [label, mat2str(val)];
    end
  elseif ischar(val) || isstring(val)
    line = [label, char(string(val))];
  else
    line = [label, sprintf('<%s>', class(val))];
  end
end

function emit(fid, fmt, varargin)
  if nargin < 2
    line = fmt;
  else
    line = sprintf(fmt, varargin{:});
  end
  fprintf('%s\n', line);
  fprintf(fid, '%s\n', line);
end

function ib = local_band_indices(E, nb)
  ib = (1:nb)';
  if isfield(E, 'config_isdf') && isstruct(E.config_isdf) ...
      && isfield(E.config_isdf, 'SYSTEM')
    sys = E.config_isdf.SYSTEM;
    if isfield(sys, 'energy_band_index_min')
      ib = double(sys.energy_band_index_min) + (0:nb - 1)';
    end
  end
end

function tag = local_case_tag(case_dir)
  [~, tag] = fileparts(case_dir);
  if isempty(tag)
    tag = 'ex_benchmark';
  end
end

function E = local_plot_ex_benchmark(E, case_dir, tag)
  E.plot_png = '';
  E.plot_fig = '';
  E.log_dir = fullfile(case_dir, 'log');

  if ~isfield(E, 'Ex_dir') || isempty(E.Ex_dir)
    warning('report_ex_benchmark:plot', '[%s] Ex empty; skip visualization.', tag);
    return
  end

  Ex_isdf = E.Ex_isdf;
  Ex_dir = E.Ex_dir;
  dEx = E.dEx;
  if ~isequal(size(Ex_isdf), size(Ex_dir), size(dEx))
    warning('report_ex_benchmark:plot', '[%s] Ex size mismatch; skip visualization.', tag);
    return
  end

  if ~isfolder(E.log_dir)
    mkdir(E.log_dir);
  end

  ts = datestr(now, 'yyyymmdd_HHMMSS');
  png_file = fullfile(E.log_dir, sprintf('%s_ex_benchmark_%s.png', tag, ts));
  fig_file = fullfile(E.log_dir, sprintf('%s_ex_benchmark_%s.fig', tag, ts));

  fig_h = figure('Name', sprintf('%s Ex benchmark', tag), 'Color', 'w', ...
    'Position', [120, 120, 1200, 760], 'Visible', 'off');

  subplot(2, 2, 1, 'Parent', fig_h);
  imagesc(Ex_dir);
  axis tight;
  colorbar;
  title('Ex\_dir (eV)');
  xlabel('ik');
  ylabel('ib');

  subplot(2, 2, 2, 'Parent', fig_h);
  imagesc(Ex_isdf);
  axis tight;
  colorbar;
  title('Ex\_isdf (eV)');
  xlabel('ik');
  ylabel('ib');

  subplot(2, 2, 3, 'Parent', fig_h);
  imagesc(dEx);
  axis tight;
  colorbar;
  title('dEx = Ex\_isdf - Ex\_dir (eV)');
  xlabel('ik');
  ylabel('ib');

  subplot(2, 2, 4, 'Parent', fig_h);
  histogram(dEx(:), 40);
  grid on;
  title('Histogram of dEx');
  xlabel('dEx (eV)');
  ylabel('Count');

  sgtitle(sprintf('%s: ISDF vs Dense Exchange', tag));

  try
    exportgraphics(fig_h, png_file, 'Resolution', 220);
  catch
    saveas(fig_h, png_file);
  end
  savefig(fig_h, fig_file);

  E.plot_png = png_file;
  E.plot_fig = fig_file;
  fprintf('\n[%s] Saved visualization:\n  %s\n  %s\n', tag, png_file, fig_file);

  if usejava('desktop')
    set(fig_h, 'Visible', 'on');
  else
    close(fig_h);
  end
end

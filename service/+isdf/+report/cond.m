% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/04 ZZ

function cond(action, varargin)
%COND  Numeric conditioning diagnostics to named OF (filename_map.cond_report).
%
%   isdf.report.cond('init')
%   isdf.report.cond('adaptive', coarse_id, idnew, desc, loss_before_step1, loss_after, final_nisdf)
%   isdf.report.cond('gen_tildevq', id)
%
% See +report/NAMING.md.

  action = lower(strtrim(char(string(action))));
  def = filename_map();
  of = def.cond_report;
  how = ['o ' of];
  ban = repmat('=', 1, 78);

  switch action
    case 'init'
      fpath = fullfile(def.isdf_report_dir, of);
      output.open(of, fpath, 'w');
      output.msg(how, '=== isdf.report.cond ===');
      output.msg(how, 'Generated: %s', datestr(now, 31));
      output.msg(how, 'dir: %s', def.isdf_report_dir);
      output.msg(how, '');

    case 'adaptive'
      coarse_id = varargin{1};
      idnew = varargin{2};
      desc = varargin{3};
      loss_before_step1 = varargin{4};
      loss_after = varargin{5};
      final_nisdf = varargin{6};
      output.msg(how, ban);
      output.msg(how, '--- adaptiveisdf ---');
      output.msg(how, 'origin_id: %d,     new_id: %d,     desc: %s,     final_nisdf: %d', ...
        round(double(coarse_id)), round(double(idnew)), char(string(desc)), round(double(final_nisdf)));
      output.msg(how, 'loss_before_step1: %.4e,     loss_after_complete: %.4e', ...
        double(loss_before_step1), double(loss_after));
      output.msg(how, '');

    case 'gen_tildevq'
      id = varargin{1};
      isdf_data = isdf.get(id);
      q = lattice.manager('q', 'get');
      nq = q.nibz;
      nisdf = double(isdf_data.nisdf);
      s_cut = double(isdf_data.svd_s_cut);
      facs = isdf_data.CCHq_trunc_factors;

      output.msg(how, ban);
      output.msg(how, '--- gen_tildeVq ---');
      output.msg(how, 'id: %d  desc: %s  nisdf: %d  s_cut: %.4e', ...
        round(double(id)), char(string(isdf_data.desc)), round(nisdf), s_cut);
      output.msg(how, '');

      for iqibz = 1:nq
        fac = facs{iqibz};
        Nkeep = double(fac.N_keep);
        lam = fac.Lambda_trunc(:);
        sigma_max = double(lam(1));
        sigma_min = double(lam(end));
        if Nkeep < nisdf
          output.msg(how, '  iqibz: %d  l_keep: %d  cond > %.4e (truncated)', ...
            iqibz, round(Nkeep), 1e12);
        else
          output.msg(how, '  iqibz: %d  l_keep: %d  cond: %.4e', ...
            iqibz, round(Nkeep), sigma_max / sigma_min);
        end
        output.msg(how, '    sigma_max: %.4e  sigma_min: %.4e  s_cut: %.4e', ...
          sigma_max, sigma_min, s_cut);
      end
      output.msg(how, '');

    otherwise
      output.err('Unknown action ''%s''.', action);
  end
end

% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/04 ZZ

function cond(action, varargin)
%COND  Numeric conditioning diagnostics to named OF o-ISDF_cond.
%
%   isdf.report.cond('init')
%   isdf.report.cond('adaptive', coarse_id, idnew, desc, loss_before_step1, loss_after, final_nisdf)
%   isdf.report.cond('gen_tildevq', id)
%
% File name from filename_map().cond_report. See +report/NAMING.md.

  action = lower(strtrim(char(string(action))));

  switch action
    case 'init'
      def = filename_map();
      output.open('cond_report', fullfile(pwd, def.cond_report), 'w');
      output.msg('o cond_report', '=== isdf numerical_cond_report ===');
      output.msg('o cond_report', 'Generated: %s', datestr(now, 31));
      output.msg('o cond_report', 'pwd: %s', pwd);
      output.msg('o cond_report', '');

    case 'adaptive'
      coarse_id = varargin{1};
      idnew = varargin{2};
      desc = varargin{3};
      loss_before_step1 = varargin{4};
      loss_after = varargin{5};
      final_nisdf = varargin{6};
      output.msg('o cond_report', '--- adaptiveisdf ---');
      output.msg('o cond_report', 'origin_id: %d,     new_id: %d,     desc: %s,     final_nisdf: %d', ...
        round(double(coarse_id)), round(double(idnew)), char(string(desc)), round(double(final_nisdf)));
      output.msg('o cond_report', 'loss_before_step1: %.4e,     loss_after_complete: %.4e', ...
        double(loss_before_step1), double(loss_after));
      output.msg('o cond_report', '');

    case 'gen_tildevq'
      id = varargin{1};
      isdf_data = isdf.get(id);
      q = lattice.manager('q', 'get');
      nq = q.nibz;
      nisdf = double(isdf_data.nisdf);
      s_cut = double(isdf_data.svd_s_cut);
      facs = isdf_data.CCHq_trunc_factors;

      output.msg('o cond_report', '--- gen_tildeVq ---');
      output.msg('o cond_report', 'id: %d  desc: %s  nisdf: %d  s_cut: %.4e', ...
        round(double(id)), char(string(isdf_data.desc)), round(nisdf), s_cut);
      output.msg('o cond_report', '');

      for iqibz = 1:nq
        fac = facs{iqibz};
        Nkeep = double(fac.N_keep);
        lam = fac.Lambda_trunc(:);
        sigma_max = double(lam(1));
        sigma_min = double(lam(end));
        if Nkeep < nisdf
          output.msg('o cond_report', '  iqibz: %d  l_keep: %d  cond > %.4e (truncated)', ...
            iqibz, round(Nkeep), 1e12);
        else
          output.msg('o cond_report', '  iqibz: %d  l_keep: %d  cond: %.4e', ...
            iqibz, round(Nkeep), sigma_max / sigma_min);
        end
        output.msg('o cond_report', '    sigma_max: %.4e  sigma_min: %.4e  s_cut: %.4e', ...
          sigma_max, sigma_min, s_cut);
      end
      output.msg('o cond_report', '');

    otherwise
      error('isdf:report:cond:action', 'Unknown action ''%s''.', action);
  end
end

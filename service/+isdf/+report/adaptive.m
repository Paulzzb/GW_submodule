% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/04 ZZ

function fpath = adaptive(r)
%ADAPTIVE  Write adaptive ISDF phase-1 report.
%
%   fpath = isdf.report.adaptive(report)
%
% Always under filename_map().isdf_report_dir / sprintf(adaptive_report, id).
% See +report/NAMING.md.

  def = filename_map();
  cid = double(r.coarse_isdf_id);
  of = sprintf(def.adaptive_report, cid);
  how = ['o ' of];
  fpath = fullfile(def.isdf_report_dir, of);
  output.open(of, fpath, 'w');

  output.msg(how, '=== adaptiveisdf phase-1 report (coarse ISDF id = %d) ===', cid);
  output.msg(how, 'Generated: %s', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
  output.msg(how, ['Elapsed wall time (tic/toc): from function start through rebuild, ', ...
    'orbit summary, and this report: %.9f s'], double(r.elapsed_phase1_seconds));
  output.msg(how, '');

  output.msg(how, '=== Adaptive ISDF Update Report ===');
  if isfield(r, 'adaptive_backend')
    output.msg(how, 'Backend            : %s', char(string(r.adaptive_backend)));
  end
  if isfield(r, 'adaptive_arithmetic')
    output.msg(how, 'Arithmetic         : %s', char(string(r.adaptive_arithmetic)));
  end
  if isfield(r, 'isdf_desc')
    output.msg(how, 'ISDF desc          : %s', char(string(r.isdf_desc)));
  end
  output.msg(how, 'Coarse ISDF id     : %d', cid);
  output.msg(how, 'Initial Nisdf      : %d', int32(r.initial_nisdf));
  output.msg(how, 'Initial bundle size: %d', int32(r.n_bundle_initial));
  output.msg(how, 'Added points       : %d', int32(r.added_nisdf));
  output.msg(how, 'Final Nisdf        : %d', int32(r.final_nisdf));
  output.msg(how, 'Final bundle size  : %d', int32(r.n_bundle_final));
  output.msg(how, 'Iterations         : %d', int32(r.iterations));
  output.msg(how, 'Initial loss       : %.8e', r.loss_initial);
  output.msg(how, 'Final loss         : %.8e', r.loss_final);
  output.msg(how, 'Relative loss      : %.8e', r.relative_loss);
  output.msg(how, 'Threshold          : %.8e', r.threshold);
  output.msg(how, 'num_add            : %d', int32(r.num_add));
  output.msg(how, 'candidate ratio    : %.8e', r.candidate_ratio);
  output.msg(how, 'ISDF ratio         : %.8e', r.isdf_ratio);
  output.msg(how, 'max cond           : %.8e', r.max_cond);
  if isfield(r, 'use_cond_guard')
    output.msg(how, 'use cond guard     : %d', logical(r.use_cond_guard));
  end
  output.msg(how, 'Nisdf cap          : %.8e', r.nmu_cap);
  output.msg(how, 'Naddmax            : %d', int32(r.naddmax));
  output.msg(how, 'param source       : %s', char(string(r.param_source)));
  output.msg(how, 'Stop reason        : %s', char(string(r.stop_reason)));
  if isfield(r, 'schur_skips')
    output.msg(how, 'Schur skips        : %d', int32(r.schur_skips));
  end
  output.msg(how, 'Stepwise loss/loss0 after each +Nadd=%3d update:', int32(r.num_add));
  n_it = int32(r.iterations);
  lh = r.loss_history;
  rlh = r.relative_loss_history;
  if n_it == 0
    output.msg(how, '  Step 0: loss/loss0 = %.8e', rlh(1));
  else
    for istep = 1:double(n_it)
      output.msg(how, '  Step %d: loss = %.8e, loss/loss0 = %.8e', ...
        istep, lh(istep + 1), rlh(istep + 1));
    end
  end
  output.msg(how, '===================================');
  output.msg(how, '');
  output.msg(how, '=== end adaptiveisdf phase-1 report ===');

  output.close(of);
end

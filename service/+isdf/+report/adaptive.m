% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/04 ZZ

function fpath = adaptive(r, outDir)
%ADAPTIVE  Write adaptive ISDF phase-1 report.
%
%   fpath = isdf.report.adaptive(report)
%   fpath = isdf.report.adaptive(report, outDir)   % default: pwd
%
% Disk name: filename_map().adaptive_report -> o-ISDF_adaptive_id%d
% See +report/NAMING.md.

  if nargin < 2 || isempty(outDir)
    outDir = pwd;
  else
    outDir = char(string(outDir));
  end
  if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
  end

  cid = double(r.coarse_isdf_id);
  def = filename_map();
  fpath = fullfile(outDir, sprintf(def.adaptive_report, cid));
  output.open('adaptive_report', fpath, 'w');

  output.msg('o adaptive_report', '=== adaptiveisdf phase-1 report (coarse ISDF id = %d) ===', cid);
  output.msg('o adaptive_report', 'Generated: %s', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
  output.msg('o adaptive_report', ['Elapsed wall time (tic/toc): from function start through rebuild, ', ...
    'orbit summary, and this report: %.9f s'], double(r.elapsed_phase1_seconds));
  output.msg('o adaptive_report', '');

  output.msg('o adaptive_report', '=== Adaptive ISDF Update Report ===');
  if isfield(r, 'adaptive_backend')
    output.msg('o adaptive_report', 'Backend            : %s', char(string(r.adaptive_backend)));
  end
  if isfield(r, 'adaptive_arithmetic')
    output.msg('o adaptive_report', 'Arithmetic         : %s', char(string(r.adaptive_arithmetic)));
  end
  if isfield(r, 'isdf_desc')
    output.msg('o adaptive_report', 'ISDF desc          : %s', char(string(r.isdf_desc)));
  end
  output.msg('o adaptive_report', 'Coarse ISDF id     : %d', cid);
  output.msg('o adaptive_report', 'Initial Nisdf      : %d', int32(r.initial_nisdf));
  output.msg('o adaptive_report', 'Initial bundle size: %d', int32(r.n_bundle_initial));
  output.msg('o adaptive_report', 'Added points       : %d', int32(r.added_nisdf));
  output.msg('o adaptive_report', 'Final Nisdf        : %d', int32(r.final_nisdf));
  output.msg('o adaptive_report', 'Final bundle size  : %d', int32(r.n_bundle_final));
  output.msg('o adaptive_report', 'Iterations         : %d', int32(r.iterations));
  output.msg('o adaptive_report', 'Initial loss       : %.8e', r.loss_initial);
  output.msg('o adaptive_report', 'Final loss         : %.8e', r.loss_final);
  output.msg('o adaptive_report', 'Relative loss      : %.8e', r.relative_loss);
  output.msg('o adaptive_report', 'Threshold          : %.8e', r.threshold);
  output.msg('o adaptive_report', 'num_add            : %d', int32(r.num_add));
  output.msg('o adaptive_report', 'candidate ratio    : %.8e', r.candidate_ratio);
  output.msg('o adaptive_report', 'max add frac       : %.8e', r.max_add_frac);
  output.msg('o adaptive_report', 'ISDF ratio         : %.8e', r.isdf_ratio);
  output.msg('o adaptive_report', 'max cond           : %.8e', r.max_cond);
  if isfield(r, 'use_cond_guard')
    output.msg('o adaptive_report', 'use cond guard     : %d', logical(r.use_cond_guard));
  end
  output.msg('o adaptive_report', 'Nisdf cap          : %.8e', r.nmu_cap);
  output.msg('o adaptive_report', 'Naddmax            : %d', int32(r.naddmax));
  output.msg('o adaptive_report', 'param source       : %s', char(string(r.param_source)));
  output.msg('o adaptive_report', 'Stop reason        : %s', char(string(r.stop_reason)));
  if isfield(r, 'schur_skips')
    output.msg('o adaptive_report', 'Schur skips        : %d', int32(r.schur_skips));
  end
  output.msg('o adaptive_report', 'Stepwise loss/loss0 after each +Nadd=%3d update:', int32(r.num_add));
  n_it = int32(r.iterations);
  lh = r.loss_history;
  rlh = r.relative_loss_history;
  if n_it == 0
    output.msg('o adaptive_report', '  Step 0: loss/loss0 = %.8e', rlh(1));
  else
    for istep = 1:double(n_it)
      output.msg('o adaptive_report', '  Step %d: loss = %.8e, loss/loss0 = %.8e', ...
        istep, lh(istep + 1), rlh(istep + 1));
    end
  end
  output.msg('o adaptive_report', '===================================');
  output.msg('o adaptive_report', '');
  output.msg('o adaptive_report', '=== end adaptiveisdf phase-1 report ===');

  output.close('adaptive_report');
end

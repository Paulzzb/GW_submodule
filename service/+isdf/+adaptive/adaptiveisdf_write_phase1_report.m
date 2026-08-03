% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/03

function fpath = adaptiveisdf_write_phase1_report(r)
%ADAPTIVEISDF_WRITE_PHASE1_REPORT  Write adaptive loop + orbit summary.
%
%   fpath = isdf.adaptive.adaptiveisdf_write_phase1_report(report)
%
% Writes adaptiveisdf_id<coarse_isdf_id>.txt in pwd. Shared by single/double
% adaptive backends (precision-agnostic formatting of the report struct).

  cid = double(r.coarse_isdf_id);
  fname = sprintf('adaptiveisdf_id%d.txt', cid);
  fpath = fullfile(pwd, fname);
  fid = fopen(fpath, 'w');
  if fid < 0
    error('adaptiveisdf:reportOpen', 'Cannot open for write: %s', fpath);
  end
  oc = onCleanup(@() fclose(fid));

  fprintf(fid, '=== adaptiveisdf phase-1 report (coarse ISDF id = %d) ===\n', cid);
  fprintf(fid, 'Generated: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
  fprintf(fid, ['Elapsed wall time (tic/toc): from function start through rebuild, ', ...
    'orbit summary, and this report: %.9f s\n\n'], double(r.elapsed_phase1_seconds));

  fprintf(fid, '\n=== Adaptive ISDF Update Report ===\n');
  if isfield(r, 'adaptive_backend')
    fprintf(fid, 'Backend            : %s\n', char(string(r.adaptive_backend)));
  end
  if isfield(r, 'adaptive_arithmetic')
    fprintf(fid, 'Arithmetic         : %s\n', char(string(r.adaptive_arithmetic)));
  end
  if isfield(r, 'isdf_desc')
    fprintf(fid, 'ISDF desc          : %s\n', char(string(r.isdf_desc)));
  end
  fprintf(fid, 'Coarse ISDF id     : %d\n', cid);
  fprintf(fid, 'Initial Nisdf      : %d\n', int32(r.initial_nisdf));
  fprintf(fid, 'Initial bundle size: %d\n', int32(r.n_bundle_initial));
  fprintf(fid, 'Added points       : %d\n', int32(r.added_nisdf));
  fprintf(fid, 'Final Nisdf        : %d\n', int32(r.final_nisdf));
  fprintf(fid, 'Final bundle size  : %d\n', int32(r.n_bundle_final));
  fprintf(fid, 'Iterations         : %d\n', int32(r.iterations));
  fprintf(fid, 'Initial loss       : %.8e\n', r.loss_initial);
  fprintf(fid, 'Final loss         : %.8e\n', r.loss_final);
  fprintf(fid, 'Relative loss      : %.8e\n', r.relative_loss);
  fprintf(fid, 'Threshold          : %.8e\n', r.threshold);
  fprintf(fid, 'num_add            : %d\n', int32(r.num_add));
  fprintf(fid, 'candidate ratio    : %.8e\n', r.candidate_ratio);
  fprintf(fid, 'max add frac       : %.8e\n', r.max_add_frac);
  fprintf(fid, 'ISDF ratio         : %.8e\n', r.isdf_ratio);
  fprintf(fid, 'max cond           : %.8e\n', r.max_cond);
  if isfield(r, 'use_cond_guard')
    fprintf(fid, 'use cond guard     : %d\n', logical(r.use_cond_guard));
  end
  fprintf(fid, 'Nisdf cap          : %.8e\n', r.nmu_cap);
  fprintf(fid, 'Naddmax            : %d\n', int32(r.naddmax));
  fprintf(fid, 'param source       : %s\n', char(string(r.param_source)));
  fprintf(fid, 'Stop reason        : %s\n', char(string(r.stop_reason)));
  if isfield(r, 'schur_skips')
    fprintf(fid, 'Schur skips        : %d\n', int32(r.schur_skips));
  end
  fprintf(fid, 'Stepwise loss/loss0 after each +Nadd=%3d update:\n', int32(r.num_add));
  n_it = int32(r.iterations);
  lh = r.loss_history;
  rlh = r.relative_loss_history;
  if n_it == 0
    fprintf(fid, '  Step 0: loss/loss0 = %.8e\n', rlh(1));
  else
    for istep = 1:double(n_it)
      fprintf(fid, '  Step %d: loss = %.8e, loss/loss0 = %.8e\n', istep, lh(istep + 1), rlh(istep + 1));
    end
  end
  fprintf(fid, '===================================\n\n');
  fprintf(fid, '=== end adaptiveisdf phase-1 report ===\n');
end

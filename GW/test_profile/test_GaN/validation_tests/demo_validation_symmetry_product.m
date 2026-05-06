% Demo: symmetry / pair-product validation (no new SAVE folder).
% Core routines live under GW/service/symmtest (addpath below).

this_dir = fileparts(mfilename('fullpath'));
symmtest_dir = fullfile(this_dir, '..', '..', '..', 'service', 'symmtest');
addpath(symmtest_dir);

ctx = val_init_from_stage();

opts = struct();
opts.tol = 1e-5;
opts.max_bands = min(8, double(ctx.wf_data.nb));
opts.max_pairs = 24;

report_fixed_k = val_run_fixed_k_symmetry_validation(ctx, opts);
report_pair_product = val_run_pair_k_symmetry_validation(ctx, opts);

validation_report = struct();
validation_report.generated_at = datestr(now, 'yyyy-mm-dd HH:MM:SS');
validation_report.ok = report_fixed_k.ok && report_pair_product.ok;
validation_report.fixed_k = report_fixed_k;
validation_report.pair_product = report_pair_product;

fprintf('\n=== validation summary ===\n');
fprintf('overall ok: %d\n', validation_report.ok);
fprintf('fixed-k check: ok=%d, max_err=%.3e, checks=%d\n', ...
  report_fixed_k.ok, report_fixed_k.max_rel_err, report_fixed_k.n_checks);
if isfield(report_fixed_k, 'debug_matfile') && ~isempty(report_fixed_k.debug_matfile)
  fprintf('fixed-k failure dump: %s\n', report_fixed_k.debug_matfile);
end
fprintf('pair-product check: ok=%d, max_err=%.3e, checks=%d\n', ...
  report_pair_product.ok, report_pair_product.max_rel_err, report_pair_product.n_checks);

assignin('base', 'validation_report', validation_report);
fprintf('Assigned variable: validation_report\n');

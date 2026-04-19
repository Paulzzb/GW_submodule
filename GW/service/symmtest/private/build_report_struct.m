function report = build_report_struct(name, err_list, fails, tol)
report = struct();
report.name = name;
report.tol = tol;
report.n_checks = int32(numel(err_list));
report.n_fail = int32(numel(fails));
report.max_rel_err = 0.0;
report.mean_rel_err = 0.0;
report.ok = true;
report.fail_examples = {};

if ~isempty(err_list)
  report.max_rel_err = max(err_list);
  report.mean_rel_err = mean(err_list);
end

if ~isempty(fails)
  report.ok = false;
  report.fail_examples = fails(1:min(20, numel(fails)));
end
end

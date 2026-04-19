function report = run_coarse_selection_test()
% Sandbox test for isdf.gen_indices_coarse.

  fprintf('\n=== isdftest: coarse grid selection ===\n');


  [Nmu, ind_mu] = isdf.gen_indices_coarse();

  report = struct();
  report.ok = true;
  report.Nmu = double(Nmu);
  report.n_selected = numel(ind_mu);
  report.unique_selected = numel(unique(double(ind_mu)));

  fprintf('Nmu               = %d\n', report.Nmu);
  fprintf('n_selected         = %d\n', report.n_selected);
  fprintf('unique_selected    = %d\n', report.unique_selected);

  if report.n_selected ~= report.unique_selected
    warning('isdftest: duplicate indices detected in selected subset.');
    report.ok = false;
  end

  fprintf('status            = %s\n', string(report.ok));

end

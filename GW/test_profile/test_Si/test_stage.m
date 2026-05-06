% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/23

function test_stage()
  % Load relay stage from test_relay_stage.mat and restore to all modules
  
  filePath = 'test_relay_stage.mat';
  
  if ~exist(filePath, 'file')
    error('test_stage: File not found: %s', filePath);
  end
  
  fprintf('\n=== Loading Relay Stage ===\n');
  fprintf('Loading from: %s\n', filePath);
  
  % Load stage from database
  relay.stage_from_db(filePath);
  fprintf('Stage loaded into relay cache.\n');
  
  % Restore to all modules
  fprintf('\n=== Restoring to Modules ===\n');
  report = relay.restore();
  
  % Display restore report
  if report.ok
    status_str = 'OK';
  else
    status_str = 'FAILED';
  end
  fprintf('\nRestore Status: %s\n', status_str);
  fprintf('Successfully restored: %d items\n', numel(report.saved));
  fprintf('Skipped items: %d\n', numel(report.skipped));
  fprintf('Errors: %d\n', numel(report.errors));
  
  if ~isempty(report.saved)
    fprintf('\nSuccessfully restored:\n');
    for i = 1:numel(report.saved)
      fprintf('  [OK] %s\n', report.saved{i});
    end
  end
  
  if ~isempty(report.skipped)
    fprintf('\nSkipped items:\n');
    for i = 1:numel(report.skipped)
      fprintf('  [SKIP] %s\n', report.skipped{i});
    end
  end
  
  if ~isempty(report.errors)
    fprintf('\nErrors during restore:\n');
    for i = 1:numel(report.errors)
      fprintf('  [ERROR] %s\n', report.errors{i});
    end
  end
  
  fprintf('\n=== Restore Complete ===\n\n');
end

% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/01 ZZ

%TEST_OUTPUT  Minimal self-test for service/+output (Phase 0).
%
%   From repo root in MATLAB:
%     cd tests
%     test_output
%
%   Checks: how parsing destinations, section numbering, file landing,
%   verbose filtering, push/showtag prefixes.

tests_dir = fileparts(mfilename('fullpath'));
if isempty(tests_dir)
  tests_dir = pwd;
end
gw_root = fileparts(tests_dir);

cd(gw_root);
QPstartup;
addpath(tests_dir);

% Smoke artifacts land next to this script (tests/).
outdir = tests_dir;

report_path = fullfile(outdir, 'r-test.log');
log_path = fullfile(outdir, 'l-test.log');
of_path = fullfile(outdir, 'o.qp');

output.free();
output.init('report', report_path, 'log', log_path, 'verbose', 1);

% --- section + msg destinations ---
screen = evalc([ ...
  'output.section(''*'', ''Top level'');', ...
  'output.section(''+'', ''Sub level'');', ...
  'output.msg(''rs'', ''on screen and report'');', ...
  'output.msg(''r'',  ''report only %d'', 42);', ...
  'output.msg(''l'',  ''log only'');', ...
  'output.msg(''v2s'', ''should be filtered from screen'');', ...
  'output.section(''-'');' ...
  ]);

assert(~isempty(strfind(screen, '[01] Top level')), 'top section missing on screen'); %#ok<*STREMP>
assert(isempty(strfind(screen, '[01.01] Sub level')) || output.verbose() >= 2, ...
  'sub section should be silent on screen at verbose=1');
assert(~isempty(strfind(screen, 'on screen and report')), 'rs message missing on screen');
assert(isempty(strfind(screen, 'report only')), 'r-only leaked to screen');
assert(isempty(strfind(screen, 'log only')), 'l-only leaked to screen');
assert(isempty(strfind(screen, 'should be filtered')), 'v2 not filtered');

output.free();  % flush files
assert(isfile(report_path), 'report file not created');
assert(isfile(log_path), 'log file not created');

rep_txt = fileread(report_path);
log_txt = fileread(log_path);
assert(~isempty(strfind(rep_txt, 'report only 42')), 'report content missing');
assert(~isempty(strfind(rep_txt, '[01.01] Sub level')), 'sub section missing in report');
assert(~isempty(strfind(log_txt, 'log only')), 'log content missing');
assert(isempty(strfind(log_txt, 'report only')), 'r-only should not be in log');

% --- named output file ---
output.init('verbose', 1);
output.open('qp', of_path);
output.msg('o qp', '%4d %12.6f', 3, 1.5);
output.close('qp');
assert(isfile(of_path), 'named ofile missing');
of_txt = fileread(of_path);
assert(~isempty(strfind(of_txt, '1.500000')), 'named ofile content missing');

% --- push / showtag / verbose ---
output.free();
output.verbose(1);
cleanup = output.push('unit'); %#ok<NASGU>
output.showtag(true);
scr2 = evalc([ ...
  'output.msg(''v1s'', ''hello normal'');', ...
  'output.msg(''v2s'', ''hello debug'');' ...
  ]);
assert(~isempty(strfind(scr2, '[unit] hello normal')), 'push/showtag prefix missing');
assert(isempty(strfind(scr2, 'hello debug')), 'v2 should be filtered at verbose=1');
clear cleanup

output.verbose(2);
scr3 = evalc('output.msg(''v2s'', ''now debug'');');
assert(~isempty(strfind(scr3, 'now debug')), 'verbose=2 failed');

% --- warn ---
output.free();
wscr = evalc('output.warn(''check me %d'', 7);');
assert(~isempty(strfind(wscr, '[WARN] check me 7')), 'warn missing');

output.free();

fprintf('test_output: artifacts in %s\n', outdir);
fprintf('test_output: PASS\n');

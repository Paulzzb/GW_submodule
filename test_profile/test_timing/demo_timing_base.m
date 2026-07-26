% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08
%
% Demo: timing base layout (Yambo mod_TIMING / mod_TIMING_live / mod_TIMING_logo)
% without reading GW config. Uses R2008a-friendly APIs only in +timing.
%
% After structure checks, runs two tiny CPU benchmarks and stores results in
% global_list clocks (same fields as Yambo TYPE(clock): total_time, call_number).
%
% Run from repo root after QPstartup, or:
%   cd(fileparts(mfilename('fullpath')));
%   cd ../../; QPstartup; cd test_profile/test_timing;

here = fileparts(mfilename('fullpath'));
cd(here);
cd ../../;
QPstartup;
cd(fullfile('test_profile', 'test_timing'));

timing.free();
timing.driver();

tm = timing.get();
assert(isa(tm, 'timing.base.timing_m'));
assert(tm.alloc);
assert(tm.global_list.alloc);
assert(strcmp(tm.global_list.name, 'global'));
assert(tm.global_list.nclock_max == double(timing.base.timing_m.NCLOCKX));
assert(tm.global_list.nclock == 0);
assert(tm.internal_list.alloc);
assert(tm.internal_list.nclock == 1);
assert(tm.internal_list.clocks(1).alloc);
assert(strcmp(tm.internal_list.clocks(1).name, 'internal'));
assert(tm.logo.ID_logo == -1);

fprintf('\n=== timing demo: register clocks and run short work (cputime) ===\n');

% --- Clock 1: dense matmul -------------------------------------------------
list = tm.global_list;
list = list.allocate_next_clock('demo_matmul');
idx1 = double(list.nclock);
t0 = cputime;
n = 160;
for k = 1:15
  A = rand(n, n);
  C = A * A;
end
dt1 = cputime - t0;

c1 = list.clocks(idx1);
c1.call_number = int32(1);
c1.start = t0;
c1.stop = t0 + dt1;
c1.total_time = dt1;
c1.running = false;
list.clocks(idx1) = c1;
tm.global_list = list;

% --- Clock 2: FFT batch ----------------------------------------------------
list = tm.global_list;
list = list.allocate_next_clock('demo_fft');
idx2 = double(list.nclock);
t0 = cputime;
x = rand(1, 65536);
for k = 1:200
  y = fft(x);
end
dt2 = cputime - t0;

c2 = list.clocks(idx2);
c2.call_number = int32(1);
c2.start = t0;
c2.stop = t0 + dt2;
c2.total_time = dt2;
c2.running = false;
list.clocks(idx2) = c2;
tm.global_list = list;

timing.save2mod(tm);

tm = timing.get();
assert(tm.global_list.nclock == 2);
assert(strcmp(tm.global_list.clocks(1).name, 'demo_matmul'));
assert(strcmp(tm.global_list.clocks(2).name, 'demo_fft'));
assert(tm.global_list.clocks(1).total_time >= 0);
assert(tm.global_list.clocks(2).total_time >= 0);

fprintf('  %-14s  calls=%3d  total_cpu=%8.4f s\n', ...
  char(tm.global_list.clocks(1).name), double(tm.global_list.clocks(1).call_number), tm.global_list.clocks(1).total_time);
fprintf('  %-14s  calls=%3d  total_cpu=%8.4f s\n', ...
  char(tm.global_list.clocks(2).name), double(tm.global_list.clocks(2).call_number), tm.global_list.clocks(2).total_time);
fprintf('  (internal) %-10s  calls=%3d  running=%d\n', ...
  char(tm.internal_list.clocks(1).name), double(tm.internal_list.clocks(1).call_number), tm.internal_list.clocks(1).running);

fprintf('=== timing.report (console, three-rule table) ===\n');
timing.report();

log_path = fullfile(here, '_demo_timing_report.log');
fid = fopen(log_path, 'w');
if fid > 0
  fclose(fid);
end
timing.report('logfile', log_path);
fidl = fopen(log_path, 'r');
if fidl > 0
  L = fgetl(fidl);
  fclose(fidl);
  assert(ischar(L) && ~isempty(strfind(L, '----')));
end
delete(log_path);

fprintf('=== relay snapshot includes updated clocks ===\n');
stage = relay.collect();
assert(isfield(stage, 'timing'));
assert(isfield(stage.timing, 'main'));
tm_snap = stage.timing.main;
assert(tm_snap.global_list.nclock == 2);

timing.free();
timing.save2mod(tm_snap);
assert(timing.get().alloc);
assert(timing.get().global_list.nclock == 2);

timing.free();
report = service_reset_persistent();
assert(any(strcmp(report.cleared, 'timing.free')));

fprintf('demo_timing_base: OK\n');

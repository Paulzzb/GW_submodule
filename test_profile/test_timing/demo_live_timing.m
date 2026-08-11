% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08
%
% Demo: Yambo-style LIVE_timing progress bar (hash bar + elapsed/estimated-total strings).
% isdf_coarse_validate_energies uses the same timing.LIVE_timing API (nhash=20 there).
%
% From GW/test_profile/test_timing: cd ../../; QPstartup; cd test_profile/test_timing;

here = fileparts(mfilename('fullpath'));
cd(here);
cd ../../;
QPstartup;
cd(fullfile('test_profile', 'test_timing'));

timing.free();
timing.driver();

tm = timing.get();
tm.live.live_report_min_seconds = 0;
timing.save2mod(tm);

n_total = 200;
timing.LIVE_timing('demo LIVE hash bar', n_total);

for k = 1:n_total
  A = rand(30, 30);
  A * A;
  timing.LIVE_timing(1);
end

timing.LIVE_timing();

fprintf('demo_live_timing: OK\n');

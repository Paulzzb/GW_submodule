% Smoke test: driver-style adaptive backend selection (desc + threshold).
% Run from test_Si: run_adaptiveisdf_precision_smoke

here = fileparts(mfilename('fullpath'));
cd(here);
cd ../../;
QPstartup;
cd(here);

load('SAVE/config.mat', 'config');
cfg = config.ISDF;

assert(local_pick_backend('vc', 2e-6) == "single");
assert(local_pick_backend('vc', 1e-8) == "single");
assert(local_pick_backend('vn', 2e-4) == "single");
assert(local_pick_backend('vn', 1e-8) == "double");
assert(local_pick_backend('nn', 2e-4) == "single");
assert(local_pick_backend('nn', 1e-8) == "double");

fprintf('run_adaptiveisdf_precision_smoke: ALL OK\n');

function backend = local_pick_backend(isdf_type, thr)
  cutoff = 1e-6;
  switch lower(isdf_type)
    case 'vc'
      backend = "single";
    case {'vn', 'nn'}
      if thr > cutoff
        backend = "single";
      else
        backend = "double";
      end
    otherwise
      error('bad type');
  end
end

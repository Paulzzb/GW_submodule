function ctx = val_init_from_stage()
% Initialize runtime managers from staged snapshot in test_Si.

cpath = fileparts(mfilename('fullpath'));
parent_dir = fileparts(cpath);

cd ../../../
QPstartup
cd(cpath)

service_reset_persistent;
packages_reset_persistent;

addpath(parent_dir);
cleanup_obj = onCleanup(@() rmpath(parent_dir)); %#ok<NASGU>
test_stage;

ctx = struct();
ctx.symm_data = symmetry.manager('get');
ctx.k_data = lattice.manager('k', 'get');
ctx.fft_data = FFT.manager('get');
ctx.wf_data = wave_functions.manager('get');
end

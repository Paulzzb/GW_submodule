function E = run_test_profile_case(case_dir, profileflag)
%RUN_TEST_PROFILE_CASE Shared launcher for nok_si8 multik flow.
%   E = run_test_profile_case(CASE_DIR, PROFILEFLAG) prepares service/MEX
%   dependencies, then runs the same GW multik flow as
%   test_nok_si8/test_gw_cohsex_multi_k.m in CASE_DIR.

if nargin < 1 || isempty(case_dir)
    case_dir = pwd;
end
if nargin < 2 || isempty(profileflag)
    profileflag = false;
end

if ~isfolder(case_dir)
    error('run_test_profile_case:MissingCaseDir', ...
        'Case directory not found: %s', case_dir);
end

test_profile_dir = fileparts(mfilename('fullpath'));
gw_root_dir = fileparts(test_profile_dir);
service_dir = fullfile(gw_root_dir, 'service');

if ~isfolder(service_dir)
    error('run_test_profile_case:MissingServiceDir', ...
        'Service directory not found: %s', service_dir);
end

cd(service_dir);
addpath(genpath(service_dir));
rehash;

ensure_mex_kernel('isdf.adaptive.isdf_schur_rank1_mex', ...
    @() isdf.adaptive.build_schur_rank1_mex);
ensure_mex_kernel('isdf.prod_mex', @() isdf.build_prod_mex);
ensure_mex_kernel('isdf.adaptive.isdf_schur_rank1_prod_mex', ...
    @() isdf.adaptive.build_schur_rank1_prod_mex);

fprintf('isdf_schur_rank1_mex available: %d\n', ...
    ~isempty(which('isdf.adaptive.isdf_schur_rank1_mex')));
fprintf('isdf_prod_mex available: %d\n', ~isempty(which('isdf.prod_mex')));
fprintf('isdf_schur_rank1_prod_mex available: %d\n', ...
    ~isempty(which('isdf.adaptive.isdf_schur_rank1_prod_mex')));

cd(gw_root_dir);
QPstartup;
cd(case_dir);

service_reset_persistent();
packages_reset_persistent();

if profileflag
    profile clear
    profile on
end
t_in = tic;
input_driver('./test');
load(fullfile(case_dir, 'SAVE', 'config.mat'), 'config');
wall_input = toc(t_in);
if profileflag
    profile off
    local_profsave_html(profile('info'), fullfile(case_dir, 'profile_input'));
end

if profileflag
    profile clear
    profile on
end
fprintf('\n[test_nok_si8] Running qp_cohsex(config) ...\n');
t_qp = tic;
E = qp_cohsex(config);
wall_qp = toc(t_qp);
if profileflag
    profile off
    local_profsave_html(profile('info'), fullfile(case_dir, 'profile_qp'));
end

fprintf('\n[test_nok_si8] Wall clock (tic/toc):\n');
fprintf('  input_driver + load(config): %.3f s\n', wall_input);
fprintf('  qp_cohsex (Ex + COHSEX):      %.3f s\n', wall_qp);
fprintf('  total:                        %.3f s\n', wall_input + wall_qp);

fprintf('\n[test_nok_si8] size(E.Eqp)    = [%d %d]\n', size(E.Eqp, 1), size(E.Eqp, 2));
fprintf('[test_nok_si8] size(E.Ex)     = [%d %d]\n', size(E.Ex, 1), size(E.Ex, 2));
fprintf('[test_nok_si8] size(E.Esx_x)  = [%d %d]\n', size(E.Esx_x, 1), size(E.Esx_x, 2));
fprintf('[test_nok_si8] size(E.Ecoh)   = [%d %d]\n', size(E.Ecoh, 1), size(E.Ecoh, 2));
fprintf('[test_nok_si8] ||E.Eqp||_F    = %.6e\n', norm(E.Eqp, 'fro'));
fprintf('[test_nok_si8] ||E.Esx_x||_F  = %.6e\n', norm(E.Esx_x, 'fro'));
fprintf('[test_nok_si8] ||E.Ecoh||_F   = %.6e\n', norm(E.Ecoh, 'fro'));

end

function ensure_mex_kernel(symbol_name, build_fn)
if isempty(which(symbol_name))
    fprintf('Building missing MEX: %s\n', symbol_name);
    build_fn();
end
end

function local_profsave_html(s, destDir)
if isempty(s)
    warning('run_test_profile_case:profsave', ...
        'profile(''info'') is empty; skip HTML export for %s.', destDir);
    return
end
if ~isfolder(destDir)
    mkdir(destDir);
end
try
    profsave(s, destDir);
catch ME
    try
        profsave(s);
        warning('run_test_profile_case:profsave', ...
            'profsave(s, dest) failed (%s); saved to default location instead.', ...
            ME.message);
    catch ME2
        throw(ME2);
    end
end
end

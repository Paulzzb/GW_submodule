function out_dat = run_fullfreq_isdf_vs_dense()
%RUN_FULLFREQ_ISDF_VS_DENSE  Compare full-frequency Gamma results on Si8.
%   Runs both ISDF and non-ISDF full-frequency (Gamma, double-k/q) using
%   the test_nok_si8 SAVE/data.mat + SAVE/config.mat inputs, then writes
%   sres/sint and their differences to a .dat file.

case_dir = fileparts(mfilename('fullpath'));
save_dir = fullfile(case_dir, 'SAVE');
out_dat = fullfile(case_dir, 'fullfreq_isdf_vs_dense_si8.dat');

if ~isfolder(save_dir)
  error('run_fullfreq_isdf_vs_dense:MissingSaveDir', ...
    'SAVE directory not found: %s', save_dir);
end

cfg_payload = load(fullfile(save_dir, 'config.mat'), 'config');
data_payload = load(fullfile(save_dir, 'data.mat'), 'data');
config0 = cfg_payload.config;
data = data_payload.data;

config0 = local_set_fullfreq_defaults(config0);

fprintf('[fullfreq-compare] Case: %s\n', case_dir);
fprintf('[fullfreq-compare] Output: %s\n', out_dat);

% ----------------------------- ISDF run -----------------------------
config_isdf = config0;
config_isdf.ISDF.isisdf = true;
config_isdf.ISDF.compute_vc = true;
config_isdf.ISDF.compute_vn = true;
config_isdf.ISDF.compute_nn = true;

service_reset_persistent();
packages_reset_persistent();
service_driver(data, config_isdf);
config_isdf = generate_frequency([], config_isdf);

t_isdf = tic;
sres_isdf = gw_fullfreq_cd_res_Gamma(config_isdf);
sint_isdf = gw_fullfreq_cd_int_Gamma(config_isdf);
fprintf('[fullfreq-compare] ISDF run done in %.3f s\n', toc(t_isdf));

% --------------------------- non-ISDF run ---------------------------
config_dense = config0;
config_dense.ISDF.isisdf = false;

service_reset_persistent();
packages_reset_persistent();
service_driver(data, config_dense);
config_dense = generate_frequency([], config_dense);

t_dense = tic;
sres_dense = gw_fullfreq_cd_res_Gamma(config_dense);
sint_dense = gw_fullfreq_cd_int_Gamma(config_dense);
fprintf('[fullfreq-compare] non-ISDF run done in %.3f s\n', toc(t_dense));

% ------------------------------ output ------------------------------
band_index = (config0.SYSTEM.energy_band_index_min:config0.SYSTEM.energy_band_index_max).';
dsres = sres_isdf - sres_dense;
dsint = sint_isdf - sint_dense;

fid = fopen(out_dat, 'w');
if fid < 0
  error('run_fullfreq_isdf_vs_dense:OpenFailed', ...
    'Cannot open output file: %s', out_dat);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '# Si8 full-frequency Gamma comparison (ISDF vs non-ISDF)\n');
fprintf(fid, '# columns (width=20):\n');
fprintf(fid, '#');
labels = { ...
  'band', ...
  'sres_isdf_re', 'sres_isdf_im', ...
  'sres_dense_re', 'sres_dense_im', ...
  'dsres_re', 'dsres_im', ...
  'sint_isdf_re', 'sint_isdf_im', ...
  'sint_dense_re', 'sint_dense_im', ...
  'dsint_re', 'dsint_im'};
for il = 1:numel(labels)
  fprintf(fid, '%20s', labels{il});
end
fprintf(fid, '\n');

for i = 1:numel(band_index)
  fprintf(fid, ['%20d', repmat('%20.4e', 1, 12), '\n'], ...
    band_index(i), ...
    real(sres_isdf(i)), imag(sres_isdf(i)), ...
    real(sres_dense(i)), imag(sres_dense(i)), ...
    real(dsres(i)), imag(dsres(i)), ...
    real(sint_isdf(i)), imag(sint_isdf(i)), ...
    real(sint_dense(i)), imag(sint_dense(i)), ...
    real(dsint(i)), imag(dsint(i)));
end

fprintf('[fullfreq-compare] Wrote %d bands to %s\n', numel(band_index), out_dat);

end

function config = local_set_fullfreq_defaults(config)
config.FREQUENCY.frequency_dependence = 2;
config.FREQUENCY.frequency_dependence_method = 2;
config.FREQUENCY.cd_residual_method = 0;
config.FREQUENCY.cd_integration_method = 0;
% Force a concrete full-frequency setup so both sres and sint are meaningful.
config.FREQUENCY.eta = 0.1;
config.FREQUENCY.frequency_low_cutoff = 2.0;
config.FREQUENCY.delta_frequency = 1.0;
config.FREQUENCY.number_imaginary_freqs = 15;
config.FREQUENCY.cd_integration_parameter = 2.0;
end

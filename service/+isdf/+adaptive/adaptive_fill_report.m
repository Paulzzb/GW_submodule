% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/03

function report = adaptive_fill_report(report, id_old, id_new, params)
%
% ADAPTIVE_FILL_REPORT  Fill report fields obtainable from stored structs / *.get.
% 
if nargin < 1 || isempty(report)
  report = struct();
end
if nargin < 4 || isempty(params)
  error('adaptive_fill_report:params', 'params struct is required.');
end

id_old = int32(id_old);
id_new = int32(id_new);
isdf_data = isdf.get(id_old);
isdf_data_new = isdf.get(id_new);

% --- Identity / sizes from ISDF objects ---
report.coarse_isdf_id = id_old;
report.isdf_desc = char(string(isdf_data.desc));

if isfield(isdf_data_new, 'N_coarse') && ~isempty(isdf_data_new.N_coarse)
  Nisdf0 = double(isdf_data_new.N_coarse);
else
  Nisdf0 = double(isdf_data.nisdf);
end
if isfield(isdf_data_new, 'N_extra') && ~isempty(isdf_data_new.N_extra)
  Nextra = double(isdf_data_new.N_extra);
else
  Nextra = double(isdf_data_new.nisdf) - Nisdf0;
end
Nisdf1 = double(isdf_data_new.nisdf);

report.initial_nisdf = int32(Nisdf0);
report.added_nisdf = int32(Nextra);
report.final_nisdf = int32(Nisdf1);

% --- Params ---
report.threshold = params.threshold;
report.num_add = int32(params.num_add);
report.candidate_ratio = params.candidate_ratio;
report.max_add_frac = params.max_add_frac;
report.isdf_ratio = params.isdf_ratio;
report.max_cond = params.max_cond;
report.use_cond_guard = logical(params.use_cond_guard);
report.param_source = params.source;

% --- Caps (params + lattice + nrange) ---
k_data = lattice.manager('k', 'get');
if ~isempty(isdf_data.nrange1) && ~isempty(isdf_data.nrange2)
  n1 = double(numel(isdf_data.nrange1));
  n2 = double(numel(isdf_data.nrange2));
else
  n1 = 0;
  n2 = 0;
end
report.nmu_cap = params.isdf_ratio * sqrt(double(k_data.nbz)) * params.max_add_frac * sqrt(n1 * n2);
report.naddmax = int32(max(0, ceil(report.nmu_cap - Nisdf0)));

% --- Bundle sizes (old vs refreshed new) ---
if isfield(isdf_data, 'bundle_struct') && isstruct(isdf_data.bundle_struct) ...
    && isfield(isdf_data.bundle_struct, 'N_bundle')
  report.n_bundle_initial = double(isdf_data.bundle_struct.N_bundle);
else
  report.n_bundle_initial = 0;
end
if isfield(isdf_data_new, 'bundle_struct') && isstruct(isdf_data_new.bundle_struct) ...
    && isfield(isdf_data_new.bundle_struct, 'N_bundle')
  report.n_bundle_final = double(isdf_data_new.bundle_struct.N_bundle);
else
  report.n_bundle_final = 0;
end
end

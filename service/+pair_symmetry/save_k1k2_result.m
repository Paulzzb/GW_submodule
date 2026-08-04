% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/19 ZZ

function save_k1k2_result(save_path, k1k2_representation, k1k2_mapping, marked, tol, nk, nsym, validation_report)
% Print concise summary and save mapping artifacts for later inspection.
%
% This helper is intentionally lightweight: it does not rebuild anything,
% it only summarizes the current result and persists raw arrays to MAT.

nrep = size(k1k2_representation, 1);
covered_ratio = nnz(marked) / numel(marked);

fprintf('\n=== k1k2 representative demo ===\n');
fprintf('nk = %d, nsym = %d\n', nk, nsym);
fprintf('nrep = %d\n', nrep);
fprintf('marked coverage = %.4f\n', covered_ratio);

if nargin >= 8 && isstruct(validation_report) && isfield(validation_report, 'ok')
	fprintf('validation ok = %d\n', validation_report.ok);
	if isfield(validation_report, 'max_shared_g0_error')
		% When validation is available, this is the largest deviation from the
		% exact shared-G0 relation over all tested target pairs.
		fprintf('max shared-G0 error = %.3e\n', validation_report.max_shared_g0_error);
	end
else
	validation_report = struct();
end

save(save_path, 'k1k2_representation', 'k1k2_mapping', 'marked', 'tol', 'validation_report');
fprintf('Saved: %s\n', save_path);
end

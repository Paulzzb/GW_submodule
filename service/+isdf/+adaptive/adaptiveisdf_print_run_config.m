% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/03

function adaptiveisdf_print_run_config(id, isdf_data, params, threshold, num_add, ratio, ...
    adaptive_backend, adaptive_arithmetic)
%ADAPTIVEISDF_PRINT_RUN_CONFIG  Screen/report start banner for adaptive ISDF.
%
%   isdf.adaptive.adaptiveisdf_print_run_config(id, isdf_data, params, ...
%     threshold, num_add, ratio, adaptive_backend, adaptive_arithmetic)

  output.msg('nrs', '[Adaptive ISDF] start  desc=%s  id=%d  backend=%s', ...
    char(string(isdf_data.desc)), int32(id), adaptive_backend);
  output.msg('rs', '  Nisdf=%d  thr=%.3e  num_add=%d  cand_ratio=%.2f', ...
    int32(isdf_data.nisdf), threshold, int32(num_add), ratio);

  output.msg('r', '=== Adaptive ISDF start ===');
  output.msg('r', 'Backend            : %s', adaptive_backend);
  output.msg('r', 'Arithmetic         : %s', adaptive_arithmetic);
  output.msg('r', 'ISDF desc          : %s', char(string(isdf_data.desc)));
  output.msg('r', 'Coarse ISDF id     : %d', int32(id));
  output.msg('r', 'Initial Nisdf      : %d', int32(isdf_data.nisdf));
  output.msg('r', 'Threshold          : %.8e', threshold);
  output.msg('r', 'num_add            : %d', int32(num_add));
  output.msg('r', 'candidate ratio    : %.4f', ratio);
  output.msg('r', 'ISDF ratio         : %.4f', params.isdf_ratio);
  output.msg('r', 'max cond           : %.4e', params.max_cond);
  output.msg('r', 'use cond guard     : %d', logical(params.use_cond_guard));
  output.msg('r', 'weight batch size  : %d', int32(params.weight_batch_size));
  output.msg('r', 'param source       : %s', params.source);
  output.msg('r', '===========================');
end

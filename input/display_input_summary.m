% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/01 ZZ

function display_input_summary(config)
%DISPLAY_INPUT_SUMMARY  Print a short screen digest; full dump to report.
%
%   Screen (<=5 lines): key CONTROL / ISDF / FREQUENCY fields.
%   Report: complete parameter dump (former fprintf body).

  % Open r-<prefix>.log next to the job (output_dir), never inside SAVE.
  prefix = 'QP';
  if isfield(config, 'CONTROL') && isfield(config.CONTROL, 'prefix') ...
      && ~isempty(config.CONTROL.prefix)
    prefix = char(string(config.CONTROL.prefix));
  end
  out_dir = '.';
  if isfield(config, 'CONTROL') && isfield(config.CONTROL, 'output_dir') ...
      && ~isempty(config.CONTROL.output_dir)
    out_dir = char(string(config.CONTROL.output_dir));
  end
  verb = 1;
  if isfield(config, 'CONTROL') && isfield(config.CONTROL, 'log_level') ...
      && ~isempty(config.CONTROL.log_level)
    verb = config.CONTROL.log_level;
  end
  report_path = fullfile(out_dir, sprintf('r-%s.log', prefix));
  output.ensure_init('report', report_path, 'verbose', verb);

  % ---- screen digest (s only; full detail goes to report) ----
  gs_type = char(string(config.CONTROL.groundstate_type));
  gs_dir = char(string(config.CONTROL.groundstate_dir));
  stor = char(string(config.CONTROL.storage_dir));
  output.msg('ns', '========== GW Input Summary ==========');
  output.msg('s',  ' Groundstate : %s  %s', gs_type, gs_dir);
  output.msg('s',  ' Storage Dir : %s', stor);
  output.msg('s',  ' ISDF / freq : enabled=%d  freq=%d', ...
    config.ISDF.isisdf,  config.FREQUENCY.frequency_dependence);
  % ensure_init stores an absolute path; show that so cwd mistakes are obvious.
  output.msg('s',  ' Report file : %s', output.get_report_path());

  % ---- full report ----
  output.msg('nr', '========== GW Input Summary ==========');
  output.msg('r',  ' Groundstate Type     : %s', gs_type);
  output.msg('r',  ' Groundstate Dir      : %s', gs_dir);
  output.msg('r',  ' Storage Dir          : %s', stor);
  output.msg('r',  ' Output Dir           : %s', char(string(config.CONTROL.output_dir)));
  output.msg('r',  ' Output File          : %s', char(string(config.CONTROL.outfile)));
  output.msg('r',  ' Prefix               : %s', char(string(config.CONTROL.prefix)));

  output.msg('r', '----------- ISDF Settings -----------');
  output.msg('r', ' ISDF Enabled         : %d', config.ISDF.isisdf);
  output.msg('r', ' ISDF Type Ratio T1/T2/T3: %.2f / %.2f / %.2f', ...
    config.ISDF.isdf_ratio_type1, config.ISDF.isdf_ratio_type2, config.ISDF.isdf_ratio_type3);
  output.msg('r', ' Adaptive Thres T1/T2/T3: %.2e / %.2e / %.2e', ...
    config.ISDF.adaptive_threshold_type1, config.ISDF.adaptive_threshold_type2, ...
    config.ISDF.adaptive_threshold_type3);
  output.msg('r', ' Adaptive NumAdd T1/T2/T3: %d / %d / %d', ...
    int32(config.ISDF.adaptive_num_add_type1), int32(config.ISDF.adaptive_num_add_type2), ...
    int32(config.ISDF.adaptive_num_add_type3));
  output.msg('r', ' Adaptive CandR T1/T2/T3: %.2f / %.2f / %.2f', ...
    config.ISDF.adaptive_candidate_ratio_type1, config.ISDF.adaptive_candidate_ratio_type2, ...
    config.ISDF.adaptive_candidate_ratio_type3);
  output.msg('r', ' Adaptive MaxCond T1/T2/T3: %.2e / %.2e / %.2e', ...
    config.ISDF.adaptive_max_cond_type1, config.ISDF.adaptive_max_cond_type2, ...
    config.ISDF.adaptive_max_cond_type3);
  if isfield(config.ISDF, 'adaptive_use_cond_guard')
    output.msg('r', ' Adaptive CondGuard    : %d', logical(config.ISDF.adaptive_use_cond_guard));
  end
  if isfield(config.ISDF, 'validate_hf')
    output.msg('r', ' ISDF validate_hf      : %d', logical(config.ISDF.validate_hf));
  end
  if isfield(config.ISDF, 'inv_strategy')
    output.msg('r', ' ISDF inv_strategy     : %s', char(string(config.ISDF.inv_strategy)));
  end
  if isfield(config.ISDF, 'inv_param')
    output.msg('r', ' ISDF inv_param        : %.2e', double(config.ISDF.inv_param));
  end
  if isfield(config.ISDF, 'inv_ratio')
    output.msg('r', ' ISDF inv_ratio        : %.4f', double(config.ISDF.inv_ratio));
  end
  if isfield(config.ISDF, 'auto_inv_param')
    output.msg('r', ' ISDF auto_inv_param   : %d', logical(config.ISDF.auto_inv_param));
  end
  if isfield(config.ISDF, 'order')
    output.msg('r', ' ISDF order (tildeVq)  : %d', int32(config.ISDF.order));
  end

  output.msg('r', '----------- CUTOFFS Settings -----------');
  output.msg('r', ' COULOMB TRUNCATION   : %3d', config.CUTOFFS.coulomb_truncation_method);
  output.msg('r', ' TRUNCATION RADIUS    : %.2f (Ry)', config.CUTOFFS.coulomb_truncation_parameter);
  output.msg('r', ' COULOMB CUTOFF RADIUS: %.2f (Ry)', config.CUTOFFS.coulomb_cutoff);
  if config.FREQUENCY.frequency_dependence == 1
    output.msg('r', ' DENSITY CUTOFF RADIUS: %.2f (Ry)', config.CUTOFFS.density_cutoff);
  end

  output.msg('r', '-------- Frequency Settings ---------');
  output.msg('r', ' Frequency Dependence : %d', config.FREQUENCY.frequency_dependence);
  if config.FREQUENCY.frequency_dependence == 2
    output.msg('r', ' Frequency Dependence : %d', config.FREQUENCY.frequency_dependence);
    output.msg('r', 'MORE INFORMATION NEEDED TO DISPLAY');
  end

  output.msg('r', '----------- COHSEX Settings -----------');
  if isfield(config.COHSEX, 'exact_ch')
    output.msg('r', ' COHSEX exact_ch      : %d', logical(config.COHSEX.exact_ch));
  end
  if isfield(config.COHSEX, 'ex_use_which_isdf')
    output.msg('r', ' COHSEX Ex ISDF slot  : %s', char(string(config.COHSEX.ex_use_which_isdf)));
  end

  if isfield(config, 'SUPERCELL')
    output.msg('r', '----------- SUPERCELL Settings -----------');
    if isfield(config.SUPERCELL, 'is_supercell')
      output.msg('r', ' SUPERCELL enabled    : %d', logical(config.SUPERCELL.is_supercell));
    end
    if isfield(config.SUPERCELL, 'use_sc_isdf')
      output.msg('r', ' SUPERCELL use SC_ISDF: %d', logical(config.SUPERCELL.use_sc_isdf));
    end
    if isfield(config.SUPERCELL, 'sc_adaptive')
      output.msg('r', ' SUPERCELL sc_adaptive : %d', logical(config.SUPERCELL.sc_adaptive));
    end
    if isfield(config.SUPERCELL, 'k1')
      output.msg('r', ' SUPERCELL k          : [%d %d %d]', ...
        int32(config.SUPERCELL.k1), int32(config.SUPERCELL.k2), int32(config.SUPERCELL.k3));
    end
    if isfield(config.SUPERCELL, 'isdf_source_dir') && ~isempty(config.SUPERCELL.isdf_source_dir)
      output.msg('r', ' SUPERCELL ISDF source: %s', char(string(config.SUPERCELL.isdf_source_dir)));
    end
  end

  if isfield(config, 'CONTROL') && isfield(config.CONTROL, 'groundstate_type') ...
      && strcmpi(char(string(config.CONTROL.groundstate_type)), 'formal')
    output.msg('r', '----------- FORMAL benchmark Settings -----------');
    output.msg('r', ' groundstate_type   : formal (synthetic wf/Eo, nsym=1)');
    if isfield(config, 'FORMAL')
      if isfield(config.FORMAL, 'wf_seed')
        output.msg('r', ' FORMAL wf_seed      : %d', round(config.FORMAL.wf_seed));
      end
      if isfield(config.FORMAL, 'eo_e0')
        output.msg('r', ' FORMAL eo_e0        : %.6f Ry', config.FORMAL.eo_e0);
      end
      if isfield(config.FORMAL, 'eo_delta')
        output.msg('r', ' FORMAL eo_delta     : %.6f Ry', config.FORMAL.eo_delta);
      end
    end
  end

  output.msg('r', '====================================');
end

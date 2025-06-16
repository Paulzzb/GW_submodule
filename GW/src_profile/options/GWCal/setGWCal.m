% function opt = setGWCal(opt, options)
function GWCal = setGWCal(data, config)


ry2ev = 13.60569253;
TOL_ZERO = 1e-12;

% Initialize options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% IO infos
outfile = config.CONTROL.outfile;
storage_dir = config.CONTROL.storage_dir;
output_dir = config.CONTROL.output_dir;

% Frequency dependent infos

% GWCal.freq_dep = options.frequency_dependence;
freq_dep = config.FREQUENCY.frequency_dependence;
freq_dep_method = config.FREQUENCY.frequency_dependence_method;
freq_low_cutoff = config.FREQUENCY.frequency_low_cutoff;
delta_frequency = config.FREQUENCY.delta_frequency;
number_imaginary_freqs = config.FREQUENCY.number_imaginary_freqs;
eta = config.FREQUENCY.eta;
cd_integration_method = config.FREQUENCY.cd_integration_method;
cd_integration_parameter = config.FREQUENCY.cd_integration_parameter;
cd_residual_method = config.FREQUENCY.cd_residual_method;

% Some default values
add_w0_freq = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GWCal = struct('frequency', config.FREQUENCY, ...
              'outfile', outfile, ...
              'storage_dir', storage_dir, ...
              'output_dir', output_dir ...
); 

end % EOF
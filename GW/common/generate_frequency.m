% function [grid_real, coeff_real_func, grid_imag, coeff_imag_func]
function config = generate_frequency(GWinfo, config)

default_Constant = constant_map();
nameConstants = fieldnames(default_Constant);
for i = 1:numel(nameConstants)
  eval(sprintf('%s = %.16f;', nameConstants{i}, default_Constant.(nameConstants{i})));
end

res_method = config.FREQUENCY.cd_residual_method;
int_method = config.FREQUENCY.cd_integration_method;

if (config.FREQUENCY.frequency_dependence_method == 2)
  
  eta = config.FREQUENCY.eta;
  freq_cutoff = config.FREQUENCY.frequency_low_cutoff*ry2ev;
  delta_freq = config.FREQUENCY.delta_frequency*ry2ev;
  nimagfreq = config.FREQUENCY.number_imaginary_freqs;
  imagparam = config.FREQUENCY.cd_integration_parameter*ry2ev;
  grid_imag = zeros(nimagfreq, 1);
  coeff_imag_func = cell(nimagfreq, 1);


  switch res_method
    case 0 % Berkeley WAY
      tmpfreq = 0.0;
      nfreq = 0;
      % Generate the real frequency grid
      while tmpfreq < freq_cutoff
        nfreq = nfreq+1;
        tmpfreq = tmpfreq + delta_freq;
      end
      nrealfreq = nfreq;
      grid_real = zeros(nfreq, 1);
      for i = 1:nfreq
        grid_real(i) = complex(0.0 + (i-1)*delta_freq, eta);
      end

      % Generate coefficient function for real frequency grid
      coeff_real_func = cell(nfreq, 1);
      for ifreq = 1:nfreq
        coeff_real_func{ifreq} = @(x) 0;
      end
      for ifreq = 1:nfreq-1
        left = real(grid_real(ifreq));
        right = real(grid_real(ifreq+1));
        coeff_real_func{ifreq+1} = @(x) coeff_real_func{ifreq+1}(x) + ...
          (x > left && x <= right) .* (x-left) ./(right-left);
        coeff_real_func{ifreq} = @(x) coeff_real_func{ifreq}(x) + ...
          (x > left && x <= right) .* (right-x) ./(right-left);
      end
    case 1
      msg = sprintf('FREQUENCY.cd_residual_method = %d is under developing', config.FREQUENCY.cd_residual_method);
      QPerror(msg);
    otherwise
      msg = sprintf('FREQUENCY.cd_residual_method = %d not supported', config.FREQUENCY.cd_residual_method);
      QPerror(msg);
  end
  
  
  switch int_method 
    case 0 % Berkeley WAY
      % grid_imag
      for ifreq = 1:nimagfreq
        zFreq = (1.0 / nimagfreq) * (ifreq-1);
        tmpFreq = -1.0 * imagparam * (zFreq / (zFreq - 1.0));
        grid_imag(ifreq) = complex(0.0, tmpFreq);
      end
      % coeff_imag_func
      for ifreq = 1:nimagfreq
        if (ifreq == 1)
          freqStart = imag(grid_imag(ifreq));
        else
          freqStart = imag((grid_imag(ifreq - 1) + grid_imag(ifreq)) * 0.5);
        end
        if (ifreq == nimagfreq)
          freqEnd = (3*imag(grid_imag(end)) - imag(grid_imag(end)))*0.5;
        else
          freqEnd = imag((grid_imag(ifreq) + grid_imag(ifreq + 1)) * 0.5);
        end
        coeff_imag_func{ifreq} = @(x) atan(freqEnd ./ x) - atan(freqStart ./ x); 
      end
    case 1 % Gauss-Legendre
      msg = sprintf('FREQUENCY.cd_integration_method = %d is under developing', config.FREQUENCY.cd_integration_method);
      QPerror(msg);
    otherwise
      msg = sprintf('FREQUENCY.cd_integration_method = %d not supported', config.FREQUENCY.cd_integration_method);
      QPerror(msg);
  end
elseif ((config.FREQUENCY.frequency_dependence_method == 0) || ...
        (config.FREQUENCY.frequency_dependence_method == 1))
  msg = sprintf('FREQUENCY.frequency_dependence_method = %d is frequency independent.\nReturning ...', config.FREQUENCY.frequency_dependence_method);
  QPlog(msg);
  return
else
  msg = sprintf('FREQUENCY.frequency_dependence_method = %d is under developing.\n', ...
   config.FREQUENCY.frequency_dependence_method);
  QPerror(msg);
end % if frequency_dependence == 2


config.freqinfo = struct();
config.freqinfo.grid_real = grid_real;
config.freqinfo.grid_imag = grid_imag;
config.freqinfo.coeff_real_func = coeff_real_func;
config.freqinfo.coeff_imag_func = coeff_imag_func;

msg = sprintf('Frequency grid generated with %d real and %d imaginary frequencies\n', ...
              nrealfreq, nimagfreq);
QPlog(msg);
msg = sprintf('Method to generate real frequencies: %d\n', res_method);
QPlog(msg);
msg = sprintf('Method to generate imaginary frequencies: %d\n', int_method);
QPlog(msg);

end % EOF
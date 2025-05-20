function [grid_real, coeff_real_func, grid_imag, coeff_imag_func] = freqgen(GWinfo, options, method_real, method_imag)
% Description
%     Generate frequency grid for GW calculation.
%     called in function gw_fullfreq*
% Input
%     GWinfo: class GWinfo, calculated from ground-state method.
%     options: class GWOptions.
%     method_real: Method to construct interpolation points.
%     method_imag: Method to choose numerical quadrature nodes. 
% Output
%     grid_real, grid_imag: Grid points on real/imaginary axis.
%     coeff_real_func: coefficient for interpolation around real axis.
%     coeff_imag_func: a length(grid_imag) cell that contains functions that
%     \int_0^{i \infty} d\omega' (x / (x^2 + \omega'^2))  W_p(\omega')
%   = \sum_k coeff_imag_func{k}(x) W_p(\omega_k)
ry2ev = 13.60569253;

if nargin <= 2
  method_real = 0;
end

if nargin <= 3
  method_imag = 1;
end

switch method_real
  case 0
    if isfield(options.GWCal, 'nFreq')
      nFreq = options.GWCal.nFreq;
    else
      error('options does not have field nFreq. Please use gwCalculation.m to do the jobs.');
    end

    if isfield(options.GWCal, 'nfreq_imag')
      nfreq_imag = options.GWCal.nfreq_imag;
    else
      error('options does not have field nFreq. Please use gwCalculation.m to do the jobs.');
    end
    nfreq_real = nFreq - nfreq_imag;

    if isfield(options.GWCal, 'dFreqGrid')
      dFreqGrid = options.GWCal.dFreqGrid;
    else
      error('options does not have field dFreqGrid. Please use gwCalculation.m to do the jobs.');
    end
    
    if isfield(options.GWCal, 'dFreqBrd')
      dFreqBrd = options.GWCal.dFreqBrd;
    else
      error('options does not have field dFreqBrd. Please use gwCalculation.m to do the jobs.');
    end

    grid_real = zeros(nfreq_real, 1);
    grid_real = dFreqGrid(1:nfreq_real) + dFreqBrd(1:nfreq_real);

    coeff_real_func = cell(nfreq_real, 1);
    for ifreq = 1:nfreq_real;
      coeff_real_func{ifreq} = @(x) 0;
    end
    % left = real(grid_real(1));
    % right = real(grid_real(2));
    % coeff_real_func{1} = @(x) (x >= left && x < right) .* (x-left) ./(right-left);
    for ifreq = 1:nfreq_real-1
      left = real(grid_real(ifreq));
      right = real(grid_real(ifreq+1));
      coeff_real_func{ifreq+1} = @(x) coeff_real_func{ifreq+1}(x) + ...
        (x > left && x <= right) .* (x-left) ./(right-left);
      coeff_real_func{ifreq} = @(x) coeff_real_func{ifreq}(x) + ...
        (x > left && x <= right) .* (right-x) ./(right-left);
      
    end
  otherwise
    error("Method for method_real = %d is not support.", method_real);
end

switch method_imag
  case 0 
    if isfield(options.GWCal, 'nFreq')
      nFreq = options.GWCal.nFreq;
    else
      error('options does not have field nFreq. Please use gwCalculation.m to do the jobs.');
    end

    if isfield(options.GWCal, 'nfreq_imag')
      nfreq_imag = options.GWCal.nfreq_imag;
    else
      error('options does not have field nFreq. Please use gwCalculation.m to do the jobs.');
    end

    if isfield(options.GWCal, 'dFreqGrid')
      dFreqGrid = options.GWCal.dFreqGrid;
    else
      error('options does not have field dFreqGrid. Please use gwCalculation.m to do the jobs.');
    end
    
    if isfield(options.GWCal, 'dFreqBrd')
      dFreqBrd = options.GWCal.dFreqBrd;
    else
      error('options does not have field dFreqBrd. Please use gwCalculation.m to do the jobs.');
    end

    grid_imag = zeros(nfreq_imag, 1);
    coeff_imag_func = cell(nfreq_imag, 1);

    
    grid_imag = dFreqGrid(nfreq_real+1:nFreq) + dFreqBrd(nfreq_real+1:nFreq);
    grid_imag = grid_imag * ry2ev;
    for ifreq = 1:nfreq_imag
      if (ifreq == 1)
        freqStart = imag(grid_imag(ifreq));
      else
        freqStart = imag((grid_imag(ifreq - 1) + grid_imag(ifreq)) * 0.5);
      end
      if (ifreq == options.GWCal.nfreq_imag)
        freqEnd = imag(grid_imag(end));
      else
        freqEnd = imag((grid_imag(ifreq) + grid_imag(ifreq + 1)) * 0.5);
      end
      coeff_imag_func{ifreq} = @(x) atan(freqEnd ./ x) - atan(freqStart ./ x);
    end
    grid_imag = grid_imag / ry2ev;
  case 1
    if isfield(options.GWCal, 'nfreq_imag')
      nfreq_imag = options.GWCal.nfreq_imag;
    else
      error('options does not have field nFreq. Please use gwCalculation.m to do the jobs.');
    end
    
    if isfield(options.GWCal, 'int_tol')
      int_tol = options.GWCal.int_tol;
      [weight, grid] = GaussLegendre(nfreq_imag, int_tol);
    else
      [weight, grid] = GaussLegendre(nfreq_imag);
    end

    % mapping the grid from [-1, 1] to [0, +i\infty]
    xz = @(z) sqrt(-1) * (2 ./ (1-z) - 1) * ry2ev;
    grid_imag = xz(grid);
    % grid_imag = grid_imag / ry2ev;

    % Generate coeff_imag_func{ifreq}, where now it becomes
    %     coefffunc(w-e_i, w') = (w - e_i) / ( (w-e_i).^2 + w'.^2 )
    % So fixed w', 
    %     coeff{ifreq} = coeff(w-e_i; w') = @(x) x ./ (x.^2 + w'.^2)
    for ifreq = 1:length(grid_imag)
      coeff_imag_func{ifreq} = ...
          @(x) ry2ev * weight(ifreq) * 2 ./ (1-grid(ifreq)).^2 * x ./ (x.^2 + abs(grid_imag(ifreq)).^2);
    end
    grid_imag = grid_imag / ry2ev;
end


    

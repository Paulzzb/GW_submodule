% function [GWOptions, options] = setoptions(mol, options)
function opt = setGWCal(opt, mol, options)


% Get some info from options, transform the readable varibles to unreadable one

% opt.GWCal;


ry2ev = 13.60569253;
TOL_ZERO = 1e-12;

% Initialize options

opt.GWCal = [];

if ~isfield(options, 'fileName')
  opt.GWCal.fileName = 'GW_output.mat';
else
  opt.GWCal.fileName = options.fileName;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Frequency dependent infos

opt.GWCal.freq_dep = options.frequency_dependence;

if (options.frequency_dependence ~= 2)
  return;
end

if (options.frequency_dependence == 2)
  % if ~isfield(options, 'frequency_dependence_method'
%   if isfield(options, 'delta_frequency')
%     opt.GWCal.dDeltaFreq = options.delta_frequency;
%   end
%   
%   if isfield(options, 'freq_low_cutoff')
%     opt.GWCal.dFreqCutoff1 = options.freq_low_cutoff;
%   end
%   
%   if isfield(options, 'number_imaginary_freqs')
%     opt.GWCal.nfreq_imag = options.number_imaginary_freqs;
%   end
%   
%   if ~isfield(options, 'number_imagniary_freq')
%     opt.GWCal.fnfreq_imag = 15;
%   end

% The following are from Epsilon/inread.f90  
  opt.GWCal.nq = 0;
  opt.GWCal.nq0 = 0;
  opt.GWCal.nq1 = 0;
  opt.GWCal.non_uniform = false;
  opt.GWCal.freq_dep = 0;
  opt.GWCal.freq_dep_method = 2;
  opt.GWCal.nFreq = 1;
  opt.GWCal.nfreq_imag = 15;
  opt.GWCal.dInitFreq = 0.0;
  opt.GWCal.dDeltaFreq = -1.0;
  opt.GWCal.dFreqStepIncrease = 1.0;
  opt.GWCal.dFreqCutoff1 = -1.0;
  opt.GWCal.dFreqCutoff2 = -1.0;
  opt.GWCal.dBrdning = -1.0;
  opt.GWCal.stop_after_qpt = -1;
  add_w0_freq = false;
  opt.GWCal.nfreq_group = 1;
%  opt.GWCal.nSFreq = 1;
%  opt.GWCal.dInitSFreq = 0.0;
%  opt.GWCal.dDeltaSFreq = -1.0;
%  opt.GWCal.dSFreqStepIncrease = 0.0;
%  opt.GWCal.dSFreqCutoff1 = -1.0;
%  opt.GWCal.dSFreqCutoff2 = -1.0;

%  opt.GWCal.fullConvLog=0
%  opt.GWCal.icutv=TRUNC_NONE
%  opt.GWCal.iwritecoul=0
%  opt.GWCal.truncval(:)=0.d0
%  opt.GWCal.ecuts=0.0d0
%  opt.GWCal.valueq0=0
%  opt.GWCal.iqexactlyzero=0
%  opt.GWCal.ncrit=0
%  opt.GWCal.efermi_input=0.0d0
%  opt.GWCal.rfermi=true
%  opt.GWCal.gcomm=-1
%  opt.GWCal.os_opt_ffts=0
%  opt.GWCal.min_fftgrid=true
%  opt.GWCal.lin_denominator=0d0
%  opt.GWCal.nfreq_group=1
%  opt.GWCal.eqp_corrections=false
%  opt.GWCal.intraband_flag=0
%  opt.GWCal.intraband_overlap_min=0.5d0
%  opt.GWCal.protection_window=0
%  opt.GWCal.num_cond_bands_ignore=0
%  opt.GWCal.patched_sampling=false
%  opt.GWCal.qgrid = 0
%  opt.GWCal.imaginary_frequency=2.0d0*ryd
 plasmaFreq = 2.0 * ry2ev;
%  opt.GWCal.do_rpa = false
% variables for nonblocking scheme


  if isfield(options, 'verbosity')
    peinf.verbosity = options.verbosity;
  end

  if isfield(options, 'frequency_dependence')
    opt.GWCal.freq_dep = options.frequency_dependence;
  end
  
  if isfield(options, 'frequency_dependence_method')
    opt.GWCal.freq_dep_method = options.frequency_dependence_method;
  end
  
  if isfield(options, 'init_frequency')
    opt.GWCal.dInitFreq = options.init_frequency;
  end
  
  if isfield(options, 'delta_frequency')
    opt.GWCal.dDeltaFreq = options.delta_frequency;
  end
  
  if isfield(options, 'delta_frequency_step')
    opt.GWCal.dFreqStepIncrease = options.delta_frequency_step;
  end
  
  if isfield(options, 'frequency_low_cutoff')
    opt.GWCal.dFreqCutoff1 = options.frequency_low_cutoff;
  end
  
  if isfield(options, 'frequency_high_cutoff')
    opt.GWCal.dFreqCutoff2 = options.frequency_high_cutoff;
  end
  
  if isfield(options, 'number_imaginary_freqs')
    opt.GWCal.nfreq_imag = options.number_imaginary_freqs;
  end
  
  if isfield(options, 'plasma_freq')
    plasmaFreq = options.plasma_freq;
  end
  
  if isfield(options, 'broadening')
    opt.GWCal.dBrdning = options.broadening;
  end
  
%  if isfield(options, 'init_sfrequency')
%    opt.GWCal.dInitSFreq = options.init_sfrequency;
%  end
%  
%  if isfield(options, 'delta_sfrequency')
%    opt.GWCal.dDeltaSFreq = options.delta_sfrequency;
%  end
%  
%  if isfield(options, 'delta_sfrequency_step')
%    opt.GWCal.dSFreqStepIncrease = options.delta_sfrequency_step;
%  end
%  
%  if isfield(options, 'sfrequency_low_cutoff')
%    opt.GWCal.dSFreqCutoff1 = options.sfrequency_low_cutoff;
%  end
%  
%  if isfield(options, 'sfrequency_high_cutoff')
%    opt.GWCal.dSFreqCutoff2 = options.sfrequency_high_cutoff;
%  end
  
%  if isfield(options, 'number_qpoints')
%    opt.GWCal.nq = options.number_qpoints; % FHJ: deprecated
%  end
  
  if isfield(options, 'qgrid')
    opt.GWCal.qgrid = options.qgrid;
  end
  
  if isfield(options, 'skip_epsilon')
    opt.GWCal.skip_epsilon =true
  end
  
%  if isfield(options, 'skip_chi')
%    opt.GWCal.skip_chi =true
%  end
 
  if isfield(options, 'imaginary_frequency')
    opt.GWCal.imaginary_frequency = options.imaginary_frequency;
  end
  
  if isfield(options, 'add_w0_freq')
    add_w0_freq = true;
  end
  
  if isfield(options, 'nfreq_group')
    opt.GWCal.nfreq_group = options.nfreq_group;
  end
end


if (options.frequency_dependence == 2)
    % Default settings for full-frequency calculations
    if (opt.GWCal.freq_dep_method == 3)
        if (opt.GWCal.dBrdning < 0)
            opt.GWCal.dBrdning = 0.25;
        end
        if (opt.GWCal.dFreqCutoff1 < 0)
            opt.GWCal.dFreqCutoff1 = 10;
        end
    else
        if (opt.GWCal.dBrdning < 0)
            opt.GWCal.dBrdning = 2.0;
        end
        if (opt.GWCal.dFreqCutoff1 < 0)
            opt.GWCal.dFreqCutoff1 = 200;
        end
        if (opt.GWCal.dFreqCutoff2 < 0)
            opt.GWCal.dFreqCutoff2 = 4 * opt.GWCal.dFreqCutoff1;
        end
    end
    if (opt.GWCal.dDeltaFreq < 0)
        opt.GWCal.dDeltaFreq = opt.GWCal.dBrdning;
    end
%   if (opt.GWCal.freq_dep_method == 1)
%       if (opt.GWCal.dDeltaSFreq < 0)
%           opt.GWCal.dDeltaSFreq = opt.GWCal.dDeltaFreq;
%       end
%       if (opt.GWCal.dSFreqCutoff1 < 0)
%           opt.GWCal.dSFreqCutoff1 = opt.GWCal.dFreqCutoff1;
%       end
%       if (opt.GWCal.dSFreqCutoff2 < 0)
%           opt.GWCal.dSFreqCutoff2 = opt.GWCal.dSFreqCutoff1;
%       end
%   end
    if (opt.GWCal.freq_dep_method ~= 2)
        opt.GWCal.nfreq_imag = 0;
    end
end

if (abs(opt.GWCal.dDeltaFreq) > TOL_ZERO && opt.GWCal.freq_dep == 2)
    % JRD: Only use low_freq_cutoff for contour deformation
    if (opt.GWCal.freq_dep_method == 2)
        opt.GWCal.dFreqCutoff2 = opt.GWCal.dFreqCutoff1;
    end

    tmpFreq = opt.GWCal.dInitFreq;
    iFreqCounter = 0;
    if (add_w0_freq)
        iFreqCounter = iFreqCounter + 1;
    end
    freqStep = opt.GWCal.dDeltaFreq;
    while (tmpFreq <= opt.GWCal.dFreqCutoff2)
        iFreqCounter = iFreqCounter + 1;
        if (tmpFreq < opt.GWCal.dFreqCutoff1)
            tmpFreq = tmpFreq + opt.GWCal.dDeltaFreq;
        else
            freqStep = freqStep + opt.GWCal.dFreqStepIncrease;
            tmpFreq = tmpFreq + freqStep;
        end
    end

    if (opt.GWCal.freq_dep_method == 2)
        iFreqCounter = iFreqCounter + opt.GWCal.nfreq_imag;
    end

    opt.GWCal.nFreq = iFreqCounter;

    % This condition plays nicely with the condition above when nfreq_group > npes, i.e. the conditions ensure
    % the number of processors will always be re-set to a legitimate/sensible number
    % Parrallel implement, not used now.
%    if (opt.GWCal.nfreq_group > opt.GWCal.nFreq)
%        fprintf('WARNING: Number of frequency groups cannot exceed number of frequencies computed\n');
%        fprintf('Resetting nfreq_group %d to number of frequencies computed %d\n', opt.GWCal.nfreq_group, opt.GWCal.nFreq);
%        opt.GWCal.nfreq_group = peinf.npes;
%    end

    opt.GWCal.dFreqGrid = zeros(1, opt.GWCal.nFreq);
    opt.GWCal.dFreqBrd = zeros(1, opt.GWCal.nFreq);

    iFreqCounter = 0;
    if (add_w0_freq)
        iFreqCounter = iFreqCounter + 1;
        opt.GWCal.dFreqGrid(iFreqCounter) = 0.0;
        opt.GWCal.dFreqBrd(iFreqCounter) = 1.0e-4 * (0.0 + 1.0i);
    end
    tmpFreq = opt.GWCal.dInitFreq;
    freqStep = opt.GWCal.dDeltaFreq;
    while (tmpFreq <= opt.GWCal.dFreqCutoff2)
        iFreqCounter = iFreqCounter + 1;
        opt.GWCal.dFreqGrid(iFreqCounter) = tmpFreq / ry2ev;
        opt.GWCal.dFreqBrd(iFreqCounter) = opt.GWCal.dBrdning * (0.0 + 1.0i) / ry2ev;

        if (tmpFreq < opt.GWCal.dFreqCutoff1)
            tmpFreq = tmpFreq + opt.GWCal.dDeltaFreq;
        else
            freqStep = freqStep + opt.GWCal.dFreqStepIncrease;
            tmpFreq = tmpFreq + freqStep;
        end
    end

    if (opt.GWCal.freq_dep_method == 2)
        for i = 1:opt.GWCal.nfreq_imag
            iFreqCounter = iFreqCounter + 1;
            zFreq = (1.0 / opt.GWCal.nfreq_imag) * double(i - 1);
            tmpFreq = -1.0 * plasmaFreq * (zFreq / (zFreq - 1.0));
            opt.GWCal.dFreqGrid(iFreqCounter) = 0.0;
            opt.GWCal.dFreqBrd(iFreqCounter) = tmpFreq * (0.0 + 1.0i);
        end

%        if (opt.GWCal.do_rpa)
%            % MDB these are the parameters for the CC grid
%            opt.GWCal.rpa_freq_grid = zeros(1, opt.GWCal.nfreq_imag);
%            omega_CC = zeros(1, opt.GWCal.nfreq_imag);
%            opt.GWCal.rpa_freq_grid = 0.0;
%            omega_CC = 0.0;
%            for i = 1:opt.GWCal.nfreq_imag
%                omega_CC(i) = i * pi * 0.5 / opt.GWCal.nfreq_imag;
%            end
%
%            for i = 1:(opt.GWCal.nfreq_imag - 1)
%                opt.GWCal.rpa_freq_grid(i) = pi / (opt.GWCal.nfreq_imag * (sin(omega_CC(i))^2));
%            end
%            opt.GWCal.rpa_freq_grid(opt.GWCal.nfreq_imag) = pi * 0.5 / (opt.GWCal.nfreq_imag * (sin(omega_CC(opt.GWCal.nfreq_imag))^2));
%
%            for i = 1:opt.GWCal.nfreq_imag
%                omega_CC(i) = 1.0 / tan(omega_CC(i));
%            end
%            omega_CC = omega_CC * ryd;
%
%            % Copy grid
%            for i = 1:opt.GWCal.nfreq_imag
%                opt.GWCal.dFreqBrd(i + (opt.GWCal.Nfreq - opt.GWCal.nfreq_imag)) = omega_CC(opt.GWCal.nfreq_imag - i + 1) * (0.0 + 1.0i);
%            end
%            % CC grid
%        end
    end

%    if (opt.GWCal.freq_dep_method == 1)
%        tmpFreq = opt.GWCal.dInitSFreq;
%        iFreqCounter = 0;
%        freqStep = opt.GWCal.dDeltaSFreq;
%        while (tmpFreq <= opt.GWCal.dSFreqCutoff2)
%            iFreqCounter = iFreqCounter + 1;
%            if (tmpFreq < opt.GWCal.dSFreqCutoff1)
%                tmpFreq = tmpFreq + opt.GWCal.dDeltaSFreq;
%            else
%                freqStep = freqStep + opt.GWCal.dSFreqStepIncrease;
%                tmpFreq = tmpFreq + freqStep;
%            end
%        end
%
%        opt.GWCal.nSFreq = iFreqCounter;
%
%        opt.GWCal.dSFreqGrid = zeros(1, opt.GWCal.nSFreq);
%
%        tmpFreq = opt.GWCal.dInitSFreq;
%        iFreqCounter = 0;
%        freqStep = opt.GWCal.dDeltaSFreq;
%        while (tmpFreq <= opt.GWCal.dSFreqCutoff2)
%            iFreqCounter = iFreqCounter + 1;
%            opt.GWCal.dSFreqGrid(iFreqCounter) = tmpFreq;
%
%            if (tmpFreq < opt.GWCal.dSFreqCutoff1)
%                tmpFreq = tmpFreq + opt.GWCal.dDeltaSFreq;
%            else
%                freqStep = freqStep + opt.GWCal.dSFreqStepIncrease;
%                tmpFreq = tmpFreq + freqStep;
%            end
%        end
%    end
elseif (opt.GWCal.freq_dep == 2 && abs(opt.GWCal.dDeltaFreq) <= TOL_ZERO)
    error('Illegal value for Delta Frequency in full frequency calculation');
elseif (opt.GWCal.freq_dep == 3)
    opt.GWCal.nFreq = 2;
    opt.GWCal.dFreqGrid = zeros(1, opt.GWCal.nFreq);
    opt.GWCal.dFreqBrd = zeros(1, opt.GWCal.nFreq);
    opt.GWCal.dFreqGrid = 0.0;
    opt.GWCal.dFreqBrd(1) = complex(0.0, 0.0); % In Godby-Needs, the first frequency is always zero
    opt.GWCal.dFreqBrd(2) = complex(0.0, 1.0)
else
    opt.GWCal.nFreq = 1;
    opt.GWCal.dFreqGrid = zeros(1, opt.GWCal.nFreq);
    opt.GWCal.dFreqBrd = zeros(1, opt.GWCal.nFreq);
    opt.GWCal.dFreqGrid = 0.0;
    opt.GWCal.dFreqBrd = 0.0;
end

if (opt.GWCal.freq_dep == 3)
    opt.GWCal.nfreq_imag = 1;
elseif (opt.GWCal.freq_dep ~= 2)
    opt.GWCal.nfreq_imag = 0;
end
























% The following are from Sigma/inread.f90

if (options.frequency_dependence == 2)
  opt.GWCal.freqevalmin = 0;
  opt.GWCal.freqevalstep = 0.2;  
  opt.GWCal.nfreq_group = 1;
  opt.GWCal.freq_grid_shift = 2;  
  max_freq_eval = 2.0;
  opt.GWCal.cd_int_method = 0;
  if isfield(options, 'cd_integral_method')
    opt.GWCal.cd_int_method = options.cd_int_method;
  end

  if isfield(options, 'max_freq_eval')
    max_freq_eval = options.max_freq_eval;
  end
  
  if isfield(options, 'delta_frequency_eval')
    opt.GWCal.freqevalstep = options.delta_frequency_eval;
  end
  
  if isfield(options, 'init_frequency_eval')
    opt.GWCal.freqevalmin = options.init_frequency_eval;
  end
  
  if isfield(options, 'number_frequency_eval')
    opt.GWCal.nfreqeval = options.number_frequency_eval;
  end
  
  if isfield(options, 'frequency_grid_shift')
    opt.GWCal.freq_grid_shift= options.frequency_grid_shift;
  end
end

opt.GWCal.nfreqeval = 2*double(int32((max_freq_eval + TOL_ZERO) / opt.GWCal.freqevalstep)) + 1;



return %main function
end % main function

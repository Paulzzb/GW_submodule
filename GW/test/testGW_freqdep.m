flag = 0;
ii = sqrt(-1);
eta = 1e-1;
if flag
  clear all;
  close all;
  dbstop if error
  addpath(genpath('/public/home/jlyang/mahh/hefeikssolv'));
  cd /public/home/jlyang/mahh/hefeikssolv
  randn('state', 0);
  rand('state', 0);
  KSSOLV_startup;
  cd test/testGW
  si2_freqdep
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % We gives a most simple choice for options, check ~/src/ISDF/ISDF.m 
  % and ~/src/GW/COmegaCstar.m (or ~/src/GW/README, functions) for details.
  options.isGW          = true; 
  % options.frequency_dependence      = 0; 
  options.frequency_dependence      = 2; 
  options.amin          = 5.0; % Currently, this options with unit Ha to consist with the main part of GW method 
  options.nv            = mol.nel/2;
  options.nc            = 4;
  options.input         = 'kssolv';
  options.isbse         = false;
  options.isisdf        = false;
  % options.isisdf        = true;
  options.iscauchy      = true;
  options.vcrank_ratio  = 8; 
  options.vsrank_ratio  = 8; 
  options.ssrank_ratio  = 8; 
  options.fileName      = 'GWoutput.mat';
  options.nv_ener       = 4;
  options.nc_ener       = 4;
  options.nv_oper       = 4;
  options.nc_oper       = 4;
  
  options.frequency_dependence_method= 2;
  options.setting_method = 'kssolv';
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % The following parts are necessary for full frequency approximation
  % All are corresponding to inputs in epsilon.inp --> Epsilon/inread.f90
  if options.frequency_dependence == 2
  % options.delta_frequency_step    = 20.0; % default value of delta frequency, Unit eV.
    options.frequency_low_cutoff = 200; %default value, Unit: eV
    options.broadening = 0.1;
    if(options.frequency_dependence_method == 2)
      options.number_imaginary_freqs = 15;   %default number of frequency imaginary part
    %  options.delta_freq = 20.0;  %default value of delta frequency, Unit: eV
      options.plasma_freq = 2;   %default value of plasma for frequency imaginary part, Unit: Ry
      options.cd_int_method = 0; %For contour deformation calculations, specify the integration method on the imaginary axis (default is 0)
    elseif(options.full_freq_method == 1)
      error('Error: Not support now!!');
    elseif(options.full_freq_method == 0)
      options.freqevalstep = 0.2; % similar to BGW
      options.freq_cutoff_high = 4.0*options.freq_cutoff;
      options.delta_frequency_step = 1.0;
    else
      error('Error: Not support now!!');
    end
    
    options.frequency_dependence_method= 2;
    options.delta_frequency = 20;
  end
  
  % The following parts are necessary for full frequency approximation
  % All are corresponding to inputs in epsilon.inp --> Epsilon/inread.f90
  if options.frequency_dependence == 2
    options.delta_frequency_eval = 0.2; 
    options.freqevalmin = 0;
    options.number_frequency_eval = 1;
    options.frequency_grid_shift = 2;
    options.max_freq_eval = 2.0;
    options.cd_integral_method = 0;   
  end
  
  
  
  % ksinfo.dFreqBrd(1) = 0.0;
  options = GWOptions(mol, options);
  ksinfor = ksinfo(mol, options.Groundstate);

end
nameConstants = fieldnames(options.Constant);
for i = 1:numel(nameConstants)
  fieldname = nameConstants{i};
  value = options.Constant.(fieldname);    
  if ~isempty(value)
    strEval = sprintf('%s = %.16f;', fieldname, value)
    eval(strEval);
  end
end

ksinfor.Z     = ksinfor.Z * sqrt(vol);
ev    = ksinfor.ev;
Z     = ksinfor.Z;
Dcoul = spdiags(ksinfor.coulG(:,4), 0, ng, ng);
% Dcoul(1, 1) = ksinfor.coulG0;


Mgvc = zeros(ng, nv_oper * nc_oper);
if ~isempty(ksinfor.aqs)
  for ind_nv = nv-nv_oper+1:nv
    ind_nv_true = ind_nv + (nv - nv_oper);
    Mgvc(:, (1:nc_oper) + (ind_nv_true-1) * (nc_oper)) ...
     = ksinfor.aqs{ind_nv}(:, nv+1:nv+nc_oper);
  end
else
  aqs = cell(nv + nc, 1);
  aqstmp = zeros(ng, nv + nc);
%  gvec = ksinfor.gvec;
  for ind_aqs = 1:nv+nc
    aqstmp = mtxel_sigma(ind_aqs, ksinfor, options.Groundstate, ...
        1:nv+nc);
    aqs{ind_aqs} = aqstmp;
  end
  ksinfor.aqs = aqs;
  for ind_nv = nv-nv_oper+1:nv
    ind_nv_true = ind_nv + (nv - nv_oper);
    Mgvc(:, (1:nc_oper) + (ind_nv_true-1) * (nc_oper)) ...
     = aqs{ind_nv}(:, nv+1:nv+nc_oper);
  end
  clear aqs scal
end  
Mgvc = conj(Mgvc);

A = Mgvc' * Dcoul * Mgvc * 4 / vol;
eden = (kron(ev(nv-nv_oper+1:nv), ones(nc_oper,1)) ... 
                - kron(ones(nv_oper, 1), ev(nv+1:nv+nc_oper))) + eta * ii;
Eden = diag(eden);
sqrtEden = diag(eden.^(0.5));
totest1 = Eden.^2 + sqrtEden * A * sqrtEden;
totest2 = Eden.^2 + A * Eden;
eigtotest1 = eig(totest1);
[~, index1] = sort(real(eigtotest1));
eigtotest2 = eig(totest2);
[~, index2] = sort(real(eigtotest2));
[~, indexeden] = sort(real(eden));

% [eigtotest1(index1), eigtotest2(index2), eden(indexeden), sort(eig(A))]

totest = A*Eden - Eden.^2;
omega = sqrt(-eig(totest));

out1 = calculateepsilon(omega(1), eta, ksinfor, options);
eigs(out1, 2, 'smallestabs')

function out = calculateepsilon(omega, eta, ksinfor, options)
ii = sqrt(-1);

nameConstants = fieldnames(options.Constant);
for i = 1:numel(nameConstants)
  fieldname = nameConstants{i};
  value = options.Constant.(fieldname);    
  if ~isempty(value)
    strEval = sprintf('%s = %.16f;', fieldname, value);
    eval(strEval);
  end
end

ksinfor.Z     = ksinfor.Z * sqrt(vol);
ev    = ksinfor.ev;
Z     = ksinfor.Z;
Dcoul = spdiags(ksinfor.coulG(:,4), 0, ng, ng);
Mgvc = zeros(ng, nv_oper * nc_oper);
if ~isempty(ksinfor.aqs)
  for ind_nv = nv-nv_oper+1:nv
    ind_nv_true = ind_nv + (nv - nv_oper);
    Mgvc(:, (1:nc_oper) + (ind_nv_true-1) * (nc_oper)) ...
     = ksinfor.aqs{ind_nv}(:, nv+1:nv+nc_oper);
  end
else
  aqs = cell(nv + nc, 1);
  aqstmp = zeros(ng, nv + nc);
%  gvec = ksinfor.gvec;
  for ind_aqs = 1:nv+nc
    aqstmp = mtxel_sigma(ind_aqs, ksinfor, options.Groundstate, ...
        1:nv+nc);
    aqs{ind_aqs} = aqstmp;
  end
  ksinfor.aqs = aqs;
  for ind_nv = nv-nv_oper+1:nv
    ind_nv_true = ind_nv + (nv - nv_oper);
    Mgvc(:, (1:nc_oper) + (ind_nv_true-1) * (nc_oper)) ...
     = aqs{ind_nv}(:, nv+1:nv+nc_oper);
  end
  clear aqs scal
end  
Mgvc = conj(Mgvc);

out = zeros(ng, ng);

Eden = (kron(ev(nv-nv_oper+1:nv), ones(nc_oper,1)) ... 
                - kron(ones(nv_oper, 1), ev(nv+1:nv+nc_oper)));
edenDRtmp = 0.5 .* ...
( 1.0 ./ (Eden - omega + ii*eta) ... 
+ 1.0 ./ (Eden + omega + ii*eta) );
out = eye(ng) - 4 * ksinfor.coulG(:,4) .* Mgvc * (edenDRtmp .* Mgvc') / vol;
end